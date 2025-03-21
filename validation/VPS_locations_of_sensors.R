library(tidyverse)
library(readxl)
library(terra)
library(MASS)
library(gganimate)

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

rov_plot <- ggplot(rov_pts) +
  geom_spatvector(aes(col = type)) +
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
  mutate(Start_Time_GMT = as.numeric(difftime(Start_Time_GMT, as_datetime("1899-12-31 00:00:00"), units = "hours"))) %>% 
  dplyr::select(Longitude = Start_Longitude, Latitude = Start_Latitude, 
                Date, Station_ID) %>% 
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
  ggtitle("Trap locations") +
  scale_color_viridis_d(end = 0.8) +
  theme_bw() +
  facet_wrap(~Date, ncol = 1) +
  coord_fixed(1) +
  theme(axis.ticks = element_blank(), axis.text = element_blank())

ggsave("plots/VPSvalidation_camera.jpg", cam_plot, width = 4.6, height = 5)
