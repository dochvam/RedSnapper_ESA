library(tidyverse)
library(readxl)
library(terra)
library(MASS)
library(gganimate)

chicken_rock_long <- -76.14321
chicken_rock_lat <- 34.62355




VPS_folder <- "Data/SnapperMvmtAbundanceStudy/VPS_Data/VPS-ChickenRock-01-Results-20240202/results/animal/"
VPS_files <- list.files(VPS_folder)

ROV_pts_reported <- 

ROV_files <- VPS_files[grepl("^8", VPS_files)]

ROV_positions <- lapply(file.path(VPS_folder, ROV_files),
                         read_csv) %>% 
  bind_rows() %>% 
  mutate(Date = as_date(Time))

ROV_pts <- vect(ROV_positions, geom = c("Longitude", "Latitude"),
                 crs = "+proj=longlat") %>% 
  project("ESRI:102003") %>% 
  as.data.frame(geom = "XY")

ROV_pts %>% 
  arrange(Time) %>% 
  ggplot() +
  # geom_point(aes(x, y, group = Transmitter, col = as.factor(Transmitter))) +
  geom_path(aes(x, y, group = Transmitter, col = as.factor(Transmitter))) +
  facet_wrap(~Date)
