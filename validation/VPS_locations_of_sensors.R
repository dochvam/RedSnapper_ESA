library(tidyverse)
library(readxl)
library(terra)
library(tidyterra)

chicken_rock_long <- -76.14321
chicken_rock_lat <- 34.62355


#### Telemetry vs. reported locations of tagged ROV ####

VPS_folder <- "Data/SnapperMvmtAbundanceStudy/VPS_Data/VPS-ChickenRock-01-Results-20240202/results/animal/"
VPS_files <- list.files(VPS_folder)

ROV_files <- VPS_files[grepl("^8", VPS_files)]

ROV_positions <- lapply(file.path(VPS_folder, ROV_files),
                        read_csv) %>% 
  bind_rows() %>% 
  mutate(Date = as_date(Time), type = "VPS fixes")



## Retrieve ROV deployment info (locs, times) as reported

rov_dat_raw <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/ROV/ROV_NC_2023_Formatted_.xlsx")

rov_dat <- bind_rows(
  rov_dat_raw %>% 
    dplyr::select(
      Site = `Site:`,
      Date = `Date:`,
      Latitude = LAT, 
      Longitude = LONG,
      length = `Transect Length1`,
      area = `Area Surveyed1`,
      hab_complexity = `Habitat Complexity1`,
      structured_pct = `Structured Habitat %Cover1`,
      count = Count1
    ) %>% 
    mutate(tID = 1),
  rov_dat_raw %>% 
    dplyr::select(
      Site = `Site:`,
      Date = `Date:`,
      Latitude = LAT, 
      Longitude = LONG,
      length = `Transect Length2`,
      area = `Area Surveyed2`,
      hab_complexity = `Habitat Complexity2`,
      structured_pct = `Structured Habitat %Cover2`,
      count = count2
    ) %>% 
    mutate(tID = 2)
) %>% 
  mutate(density = count/area,
         type = "Reported drops")

rov_pts <- vect(bind_rows(rov_dat, ROV_positions), 
                geom = c("Longitude", "Latitude"),
                crs = "+proj=longlat")

rov_plot <- ggplot() +
  geom_spatvector(data = rov_pts, aes(col = type)) +
  ggtitle("ROV locations") +
  scale_color_viridis_d(end = 0.8) +
  theme_minimal()
ggsave("plots/VPSvalidation_rov.jpg", rov_plot, width = 4.6, height = 2.5)


#### Telemetry vs. reported locations of tagged traps ####

VPS_folder <- "Data/SnapperMvmtAbundanceStudy/VPS_Data/VPS-ChickenRock-01-Results-20240202/results/animal/"
VPS_files <- list.files(VPS_folder)

cam_files <- VPS_files[grepl("^2058", VPS_files)]
cam_positions <- lapply(file.path(VPS_folder, cam_files),
                        read_csv) %>% 
  bind_rows() %>% 
  mutate(Date = as_date(Time), type = "VPS fixes")

stations <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/cameratraps/RS_2023_full_reads_all_three_cameratrap_dates.xlsx", sheet = "StationData") %>% 
  dplyr::select(Longitude = Start_Longitude, Latitude = Start_Latitude, 
                Date, Station_ID, Start_Time_GMT) %>% 
  mutate(type = "Reported drops") %>% 
  mutate(Longitude = as.numeric(Longitude)) %>% 
  distinct()

# Which camera is it?
# Assign each VPS fix to whichever camera was 
#   (A) deployed on the same day
#   (B) Nearest

cam_pos_pts <- vect(cam_positions, 
                    geom = c("Longitude", "Latitude"),
                    crs = "+proj=longlat") %>% 
  project("ESRI:102003") %>% 
  as.data.frame(geom = "XY")
stations_pts <- vect(stations, 
                    geom = c("Longitude", "Latitude"),
                    crs = "+proj=longlat") %>% 
  project("ESRI:102003") %>% 
  as.data.frame(geom = "XY")

for (i in 1:nrow(cam_pos_pts)) {
  today_drops <- stations_pts %>% filter(Date == as_date(cam_pos_pts$Date[i]))
  distances <- sqrt((cam_pos_pts$x[i]  - today_drops$x)^2 +
                    (cam_pos_pts$y[i]  - today_drops$y)^2)
    
  
  cam_pos_pts$dist[i] <- min(distances)
  cam_pos_pts$Station_ID[i] <- today_drops$Station_ID[cam_pos_pts$dist[i] == distances]
}


cam_pts <- bind_rows(stations_pts, cam_pos_pts)


cam_plot <- ggplot(cam_pts, aes(x, y, group = Station_ID)) +
  geom_point(aes(col = type)) +
  geom_line() +
  # geom_label(aes(x, y, label =Station_ID)) +
  ggtitle("Trap locations") +
  scale_color_viridis_d(end = 0.8) +
  theme_bw() +
  facet_wrap(~Date, ncol = 1) +
  coord_fixed(1) +
  theme(axis.ticks = element_blank(), axis.text = element_blank())

ggsave("plots/VPSvalidation_camera.jpg", cam_plot, width = 4.6, height = 5)




#### Create a dataset for the cameras giving their *corrected* locations ####

# For each camera, find the VPS fix closest to the recorded drop point 
# that occurred during the survey, if any
stations_pts$start_time <- difftime(stations_pts$Start_Time_GMT, as_date("1899-12-31"), units = "mins")
stations_pts$start_time[stations_pts$Date != as_date("2023-08-22")] <- stations_pts$start_time[stations_pts$Date != as_date("2023-08-22")] - 60

stations_pts$end_time <- stations_pts$start_time + 60
stations_pts$x_corrected <- NA
stations_pts$y_corrected <- NA
stations_pts$corrected <- NA

for (i in 1:nrow(stations_pts)) {
  
  vps_during_survey <- cam_pos_pts %>% 
    mutate(time_since_SOD = difftime(Time, Date, units = "mins")) %>% 
    filter(Date == as_date(stations_pts$Date[i]),
           time_since_SOD <= stations_pts$end_time[i],
           time_since_SOD >= stations_pts$start_time[i]
           ) %>% 
    mutate(distance = sqrt((stations_pts$x[i]  - x)^2 +
                           (stations_pts$y[i]  - y)^2))
  
  print(vps_during_survey$distance)
  
  min_dist_index <- which.min(vps_during_survey$distance)
  if (nrow(vps_during_survey) > 1 && vps_during_survey$distance[min_dist_index] < 30) {
    stations_pts$corrected[i] <- TRUE
    stations_pts$x_corrected[i] <- vps_during_survey$x[min_dist_index]
    stations_pts$y_corrected[i] <- vps_during_survey$y[min_dist_index]
  } else {
    stations_pts$corrected[i] <- FALSE
    stations_pts$x_corrected[i] <- stations_pts$x[i]
    stations_pts$y_corrected[i] <- stations_pts$y[i]
  }
}

# Convert back to long/lat
stations_pts2 <- stations_pts %>% 
  select(Date, Station_ID, start_time, x_corrected, y_corrected, corrected) %>% 
  vect(geom = c("x_corrected", "y_corrected"),
       crs = "ESRI:102003") %>% 
  project("+proj=longlat") %>% 
  as.data.frame(geom = "XY") %>% 
  rename(Longitude = x, Latitude = y)

nrow(distinct(stations_pts2, Date, Longitude, Latitude))

write_csv(stations_pts2, "intermediate/corrected_camera_stations.csv")


