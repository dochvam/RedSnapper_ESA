library(tidyverse)
library(readxl)
library(terra)
library(MASS)
library(gganimate)

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
  bind_rows() %>% 
  mutate(Date = as_date(Time))

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

fish_pts <- vect(fish_positions, geom = c("Longitude", "Latitude"),
                 crs = "+proj=longlat") %>% 
  project("ESRI:102003") %>% 
  as.data.frame(geom = "XY")

# Camera days
cam_days <- as_date(c("2023-08-12", "2023-09-08", "2023-10-30"))
camera_locs$Longitude <- as.numeric(camera_locs$Longitude)
camera_pts <- vect(camera_locs, geom = c("Longitude", "Latitude"),
                   crs = "+proj=longlat") %>% 
  project("ESRI:102003") %>% 
  as.data.frame(geom = "XY")


# For each fish
for (i in 1:nrow(all_IDs)) {
  # Does this fish have any detections on the camera days?
  this_fixes <- fish_pts %>%
    filter(FullId == all_IDs$FullId[i])
  
  for (d in 1:length(cam_days)) { if (any(this_fixes$Date == cam_days[d])) {
    this_fixes_bydate <- this_fixes %>% 
      filter(Date %in% (cam_days[d] - 3:0)) %>% 
      mutate(color = ifelse(Date == cam_days[d], "blue", "gray"))
    
    p <- ggplot() +
      geom_point(data = camera_pts[camera_pts$Date == cam_days[d], c("x", "y")], aes(x = x, y = y), color = "black", size = 1, alpha = 0.5) +
      geom_path(data = this_fixes_bydate, aes(x, y, group = 1, color = color), alpha = 0.4, size = 1) +  # Semi-transparent previous path
      geom_point(data = this_fixes_bydate, aes(x, y, group = 1), size = 3, color = "red") +  # Current position
      theme_minimal() +
      scale_color_identity() + 
      transition_reveal(Time) + # Animate along time
      ggtitle("Time: {frame_along}")  # Show time in title

    anim_save(animation = p, paste0("plots/camera_day_animations/", all_IDs$FullId[i], "_", cam_days[d], ".gif"),
              fps = 20, duration = 7, width = 350, height = 350, 
              start_pause = 20, end_pause = 40)
  }}
}

