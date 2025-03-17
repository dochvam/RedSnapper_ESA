library(tidyverse)
library(readxl)
library(terra)
library(MASS)

chicken_rock_long <- -76.14321
chicken_rock_lat <- 34.62355

####
#' What are the various goals we have with the telemetry data?
#'  > Estimate the home range size of fish for use in the uSCR model
#'    > For the uSCR we actually just want to understand behavior during a 20 minute pd?
#'  > Validate the ROV sampling locations/paths  
#'  > Validate the camera sampling locations
#'  
#'  I'll handle these in reverse order, starting with the simplest

#### Categorize the files ####

VPS_folder <- "Data/SnapperMvmtAbundanceStudy/VPS_Data/VPS-ChickenRock-01-Results-20240202/results/animal/"
VPS_files <- list.files(VPS_folder)

fish_files <- VPS_files[grepl("^12", VPS_files)]

#### Validating camera locations ####

camera_locs <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/cameratraps/RS_2023_full_reads_all_three_cameratrap_dates.xlsx", sheet = "StationData") %>% 
  filter(`Camera (A or C)` == "A") %>% 
  dplyr::select(Station_ID, Date, Time = Start_Time_GMT,
                Latitude = Start_Latitude, Longitude = Start_Longitude)



#### Analyzing fish paths ####

fish_fate <- read_xlsx("Data/SnapperMvmtAbundanceStudy/VPS_Data/Fate_assignments_from_discard_paper_vps 2023.xlsx")
fish_to_drop <- fish_fate$`Tag number`[fish_fate$`subsequent assignments based on velocity and depth data` == "Discard mortality"]
dead_fish_files <- paste0(fish_to_drop, ".csv")
fish_files <- fish_files[!fish_files %in% dead_fish_files]

fish_positions <- lapply(file.path(VPS_folder, fish_files),
                         read_csv) %>% 
  bind_rows()

fish_pts <- vect(fish_positions, geom = c("Longitude", "Latitude"),
                 crs = "+proj=longlat") %>% 
  project("ESRI:102003") %>% 
  as.data.frame(geom = "XY")

fish_fix_times <- fish_positions %>% 
  dplyr::select(FullId, Time) %>% 
  arrange(FullId, Time)

fish_fix_times$tdiff <- c(NA, difftime(fish_fix_times$Time[2:nrow(fish_fix_times)],
                                 fish_fix_times$Time[2:nrow(fish_fix_times) - 1],
                                 units = "mins"))
fish_fix_times$tdiff[c(1, 1 + which(fish_fix_times$FullId[2:nrow(fish_fix_times)] != fish_fix_times$FullId[2:nrow(fish_fix_times) - 1]))] <- NA

quantile(fish_fix_times$tdiff, probs = 0:10/10, na.rm = T)

fish_fix_times %>% 
  mutate(tdiff = ifelse(tdiff > 100, 100, tdiff)) %>% 
  ggplot() + 
  geom_histogram(aes(tdiff)) +
  theme_minimal()

all_IDs <- fish_pts %>% 
  dplyr::select(x, y, Time, FullId) %>% 
  group_by(FullId) %>% 
  mutate(time_abs = as.numeric(difftime(Time, min(Time), units = "mins"))) %>% 
  summarize(n = n(), total_time = max(time_abs), 
            total_time_days = total_time / (24 * 60),
            avg_time_btw_fixes = total_time / n)

bbox <- c(min(fish_pts$x),
          max(fish_pts$x),
          min(fish_pts$y),
          max(fish_pts$y))

kde_list <- list()
for (i in 1:nrow(all_IDs)) {
  this_dat <- fish_pts %>% 
    filter(FullId == all_IDs$FullId[i]) %>% 
    dplyr::select(x, y, Time)
  
  kde_list[[i]] <- MASS::kde2d(this_dat$x, this_dat$y, 
                               n = 250, lims = bbox)
}







result_list <- list()
num_points <- c()
t_interval <- 60*6
pb <- progress::progress_bar$new(total = nrow(all_IDs))

for (i in 1:nrow(all_IDs)) {
  pb$tick()
  this_dat <- fish_pts %>% 
    filter(FullId == all_IDs$FullId[i]) %>% 
    dplyr::select(x, y, Time) %>% 
    mutate(time_abs = as.numeric(difftime(Time, min(Time), units = "mins")))
  
  total_time <- max(this_dat$time_abs)
  
  # num_slices <- floor(total_time / t_interval)
  num_slices <- 10000
  
  ct <- 0
  random_start <- runif(num_slices, 0, total_time-t_interval)
  random_start_real <- (60 * random_start) + min(this_dat$Time)
  centroid_y <- centroid_x <- rep(NA, length(random_start))
  estimated_variance <- rep(NA, num_slices)
  for (s in 1:num_slices) {
    include <- this_dat$time_abs > random_start[s] & 
               this_dat$time_abs < (random_start[s] + t_interval)
    num_points <- c(num_points, sum(include))
    
    if (sum(include) > 1) {
      mtx <- this_dat[include, 1:2]
      centroid <- colMeans(mtx)
      
      centroid_x[s] <- centroid[1]
      centroid_y[s] <- centroid[2]
      
      estimated_variance[s] <-
        sum((mtx - matrix(centroid, nrow(mtx), 2, byrow=TRUE))^2) / (2 * (nrow(mtx)-1))
    }
  }
  result_list[[i]] <- data.frame(var = estimated_variance,
                                 random_start = random_start,
                                 random_start_real = random_start_real,
                                 x = centroid_x,
                                 y = centroid_y,
                                 hour = hour(random_start_real),
                                 time = hour(random_start_real) + minute(random_start_real)/60,
                                 FullId = all_IDs$FullId[i])
}

all_results <- bind_rows(result_list)

# Where are the data?
ggplot(all_results) +
  geom_bar(aes(as.factor(hour), group = is.na(var), fill = is.na(var)), position = "stack") +
  theme_minimal() +
  scale_fill_viridis_d(">2 fixes", end = 0.7)


# Movement distribution per hour of day
all_results %>% 
  filter(!is.na(var)) %>% 
  ggplot() +
  geom_boxplot(aes(hour, group = hour, log(var))) +
  theme_minimal() +
  scale_color_viridis_d(end = 0.8)

# Movement distribution per fish
all_results %>% 
  filter(!is.na(var)) %>% 
  ggplot() +
  geom_density(aes(log(var), group = FullId, col = FullId),
               show.legend = F) +
  theme_minimal() +
  scale_color_viridis_d(end = 0.8)

# Movement distribution per fish, in sigma
all_results %>% 
  filter(!is.na(var)) %>% 
  ggplot() +
  geom_density(aes(sqrt(var), group = FullId, col = FullId),
               show.legend = F) +
  theme_minimal() +
  scale_color_viridis_d(end = 0.8)

# Movement dist overall
bind_rows(result_list) %>% 
  ggplot() +
  geom_density(aes(log(var))) +
  theme_minimal()

mean(log(all_results$var))
median(log(all_results$var))
sd(log(all_results$var))

mean(sqrt(all_results$var))


#### Look into days/ times when bait was deployed ####

camera_locs$Datetime <- as_datetime(camera_locs$Date) + 
  60*60*hour(camera_locs$Time) + 60 * minute(camera_locs$Time)


# Fish movement on camera days
all_results %>% 
  mutate(camera_day = yday(random_start_real) %in% yday(camera_locs$Datetime)) %>% 
  ggplot() +
  geom_density(aes(log(var), group = camera_day, col = camera_day)) +
  theme_minimal()


camera_locs$Longitude <- as.numeric(camera_locs$Longitude)
camera_pts <- vect(camera_locs, geom = c("Longitude", "Latitude"),
                   crs = "+proj=longlat") %>% 
  project("ESRI:102003") %>% 
  as.data.frame(geom = "XY")

all_results$near_deploy <- F
all_results <- all_results %>% filter(!is.na(var))

for (i in 1:nrow(camera_pts)) {
  
  all_results$near_deploy[sqrt(
      (all_results$x - camera_pts$x[i])^2 +
      (all_results$y - camera_pts$y[i])^2) < 200 & 
    as.numeric(difftime(all_results$random_start_real, 
             camera_pts$Datetime[i], units = "mins")) < 60 * 6 &
    as.numeric(difftime(all_results$random_start_real, 
                        camera_pts$Datetime[i], units = "mins")) > 0
  ] <- TRUE
}

table(all_results$near_deploy)
all_results %>% 
  ggplot() +
  geom_density(aes(log(var), group = near_deploy, col = near_deploy)) +
  theme_minimal()


wilcox.test(all_results$var[all_results$near_deploy],
            all_results$var[!all_results$near_deploy])
