library(tidyverse)
library(readxl)
library(terra)
library(MASS)

set.seed(152801)

#### Load in camera data to retrieve target time interval ####
stations <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/cameratraps/RS_2023_full_reads_all_three_cameratrap_dates.xlsx", sheet = "StationData") %>% 
  mutate(Start_Time_GMT = as.numeric(difftime(Start_Time_GMT, as_datetime("1899-12-31 00:00:00"), units = "mins")))

dates <- unique(stations$Date)
t_intervals <- numeric(length(dates))
for (i in 1:length(dates)) {
  t_intervals[i] <- diff(range(stations$Start_Time_GMT[stations$Date == dates[i]]))
}
t_interval <- max(t_intervals)

#### Analyzing fish paths ####

VPS_folder <- "Data/SnapperMvmtAbundanceStudy/VPS_Data/VPS-ChickenRock-01-Results-20240202/results/animal/"
VPS_files <- list.files(VPS_folder)

fish_files <- VPS_files[grepl("^12", VPS_files)]

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
            avg_time_btw_fixes = total_time / n) %>% 
  filter(n > 1000)

bbox <- c(min(fish_pts$x),
          max(fish_pts$x),
          min(fish_pts$y),
          max(fish_pts$y))

all_results <- list()
pb <- progress::progress_bar$new(total = nrow(all_IDs))

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
  result_df <- data.frame(var = estimated_variance,
                          random_start = random_start,
                          random_start_real = random_start_real,
                          x = centroid_x,
                          y = centroid_y,
                          hour = hour(random_start_real),
                          time = hour(random_start_real) + minute(random_start_real)/60,
                          FullId = all_IDs$FullId[i])
}


write_csv(result_df, "pipeline_NC/NC_results/all_intervals_sigma.csv")

result_df <- read_csv("pipeline_NC/NC_results/all_intervals_sigma.csv") %>% 
  filter(!is.na(var))

log_sigma_summary <- result_df %>% 
  summarize(mean = mean(log(sqrt(var))),
            median = median(log(sqrt(var))),
            sd = sd(log(sqrt(var))),
            q25 = quantile(log(sqrt(var)), 0.25),
            q75 = quantile(log(sqrt(var)), 0.75))

qqnorm(scale(log(sqrt(result_df$var))))
abline(0,1)

write_csv(log_sigma_summary, "pipeline_NC/NC_results/log_sigma_estimate_NC.csv")
