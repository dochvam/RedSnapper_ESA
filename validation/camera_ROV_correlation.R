library(readxl)
library(tidyverse)
library(gridExtra)
library(terra)
library(tidyterra)

set.seed(5249938)

hardbottom <- rast("Data/Chicken_Rock_Map/ChickenRock_Classification.tif")

#### Load the camera data ####

camera_dat <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/cameratraps/RS_2023_full_reads_all_three_cameratrap_dates.xlsx") %>% 
  mutate(time = as.numeric(difftime(time, as_datetime("1899-12-31 00:00:00"), units = "sec"))) %>% 
  filter(is.na(comments) | grepl("tagged", comments)) # Filter out descent/retrieval frames

stations <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/cameratraps/RS_2023_full_reads_all_three_cameratrap_dates.xlsx", sheet = "StationData") %>% 
  mutate(Start_Time_GMT = as.numeric(difftime(Start_Time_GMT, as_datetime("1899-12-31 00:00:00"), units = "hours"))) %>% 
  mutate(Start_Longitude = as.numeric(Start_Longitude))

camera_counts <- camera_dat %>% 
  filter(!is.na(Station_ID)) %>% 
  group_by(Station_ID, date) %>% 
  summarize(max_count = max(total),
            total_count = sum(total)) %>% 
  left_join(stations, by = c("date" = "Date", "Station_ID")) %>%
  rename(Longitude = Start_Longitude, Latitude = Start_Latitude) %>% 
  mutate(current_dir = ifelse(grepl("toward", tolower(`Current Direction`)),
                              "Toward", ifelse(grepl("away", tolower(`Current Direction`)),
                                               "Away", "Perpendicular"))) %>% 
  mutate(Longitude = as.numeric(Longitude)) %>% 
  filter(`Camera (A or C)` == "A")

camera_pts <- vect(camera_counts, geom = c("Longitude", "Latitude"), 
                   crs = "+proj=longlat") %>% 
  project(crs(hardbottom))

#### Load the ROV data ####

rov_dat_raw <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/ROV/ROV_NC_2023_Formatted_.xlsx")

rov_dat <- bind_rows(
  rov_dat_raw %>% 
    select(
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
    select(
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
  group_by(Site, Longitude, Latitude, Date) %>% 
  summarize(count = sum(count), area = sum(area)) %>% 
  mutate(density = count/area)

rov_pts <- vect(rov_dat, geom = c("Longitude", "Latitude"), 
                crs = ("+proj=longlat")) %>% 
  project(crs(hardbottom))

#### For each camera, get the ID of the closest ROV survey ####

dist_mtx <- distance(camera_pts, rov_pts)

camera_counts$closest_ROV_dist <- apply(dist_mtx, 1, min)
camera_counts$closest_ROV_ID <- apply(dist_mtx, 1, which.min)
camera_counts$closest_ROV_density <- rov_pts$density[closest_ROV_ID]

camera_counts %>% 
  filter(closest_ROV_dist < 100) %>% 
  ggplot(aes(closest_ROV_density, total_count,  col = closest_ROV_dist)) + 
  geom_point() +
  geom_smooth(method = "lm")

camera_counts %>% 
  filter(closest_ROV_dist < 100) %>% 
  ggplot(aes(closest_ROV_density, max_count, col = closest_ROV_dist)) + 
  geom_point() + 
  geom_smooth(method = "lm")

summary(lm(total_count ~ closest_ROV_density, data = camera_counts))
summary(lm(max_count ~ closest_ROV_density, data = camera_counts))


