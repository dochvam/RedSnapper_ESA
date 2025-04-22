library(tidyverse)
library(readxl)
library(terra)
library(MASS)
library(suncalc)

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

all_results <- list()
pb <- progress::progress_bar$new(total = 40*nrow(all_IDs))

for (t in 1:12) {
  t_interval <- 20 * t * 2 - 20
  
  result_list <- list()
  num_points <- c()
  
  for (i in 1:nrow(all_IDs)) {
    pb$tick()
    this_dat <- fish_pts %>% 
      filter(FullId == all_IDs$FullId[i]) %>% 
      dplyr::select(x, y, Time) %>% 
      mutate(time_abs = as.numeric(difftime(Time, min(Time), units = "mins")))
    
    total_time <- max(this_dat$time_abs)
    
    num_slices <- floor(total_time / t_interval)
    
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
  
  all_results[[t]] <- bind_rows(result_list) %>% 
    mutate(t_interval = t_interval)
}


all_all_results <- bind_rows(all_results)
write_csv(all_all_results, "intermediate/all_intervals_sigma.csv")

all_all_results <- read_csv("intermediate/all_intervals_sigma.csv")

sigma_summary <- all_all_results %>% 
  filter(!is.na(var), t_interval < 800) %>% 
  group_by(t_interval) %>% 
  summarize(med = median(sqrt(var)), q25 = quantile(sqrt(var), 0.25),
            q75 = quantile(sqrt(var), 0.75))

# Plot sigma
sigma_summary %>% 
  ggplot() + 
  geom_vline(xintercept = 20, col = "red") +
  geom_vline(xintercept = 60*6, col = "red") +
  geom_pointrange(aes(t_interval, med, ymin = q25, ymax = q75)) +
  xlab("Time inteval (minutes)") + ylab("Sigma") +
  theme_minimal()

# Plot sigma vs. the expectation if we think of these as a convolution of many
# MVNs with 20-min intervals
var_20min_vec <- all_all_results %>% 
  filter(!is.na(var), t_interval == 20) %>% 
  .$var

# For each time interval, say var is a convolution of that number of normals.
# Do this 500x per.

sigma_summary_convln <- sigma_summary
for (i in 1:nrow(sigma_summary_convln)) {
  sigma_vec <- numeric(500)
  for (j in 1:500) {
    sigma_vec[j] <- sqrt(sum(sample(var_20min_vec, size = sigma_summary_convln$t_interval[i] / 20)))
  }
  
  sigma_summary_convln$med[i] <- median(sigma_vec)
  sigma_summary_convln$q25[i] <- quantile(sigma_vec, 0.25)
  sigma_summary_convln$q75[i] <- quantile(sigma_vec, 0.75)
}

bind_rows(sigma_summary_convln %>% mutate(type = "expected"),
          sigma_summary %>% mutate(type = "observed")) %>% 
  ggplot() + 
  geom_pointrange(aes(t_interval, med, ymin = q25, ymax = q75,
                      group = type, col = type), position = position_dodge(width = 10)) +
  xlab("Time inteval (minutes)") + ylab("Sigma") +
  theme_minimal()




# Plot variance
all_all_results %>% 
  filter(!is.na(var)) %>% 
  group_by(t_interval) %>% 
  summarize(mean = median(sqrt(var)), q25 = quantile(sqrt(var), 0.25),
            q75 = quantile(sqrt(var), 0.75)) %>% 
  ggplot() + 
  # geom_vline(xintercept = 20, col = "red") +
  # geom_vline(xintercept = 60*6, col = "red") +
  geom_pointrange(aes(sqrt(t_interval), mean, ymin = q25, ymax = q75)) +
  xlab("Time inteval (minutes)") + ylab("Sigma") +
  theme_minimal()

library(lme4)
library(cv)

# Log-normal model of variance
fit <- lmer(log((var)) ~ log(t_interval) + (1 | FullId), data = all_all_results)
summary(fit)

# Log-normal model of sd
all_all_results$log_sigma <- log(sqrt(all_all_results$var))
fit <- lmer(log_sigma ~ log(t_interval) + (1 | FullId), data = all_all_results)
summary(fit)

fit <- lm(log_sigma ~ log(t_interval), data = all_all_results)
summary(fit)

predict(fit, newdata = data.frame(t_interval = 60 * 6), se.fit = T)

