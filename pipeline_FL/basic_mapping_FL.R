library(readxl)
library(tidyverse)
library(gridExtra)
library(terra)
library(tidyterra)
library(measurements)

set.seed(5249938)

#### Load the VPS data ####

vps_results <- list.files("Data/SouthAtlantic/VPS_results/animal/",
                          full.names = TRUE) %>% 
  lapply(read_csv, progress = F, col_types = list(Transmitter = col_character())) %>% 
  bind_rows() %>% 
  dplyr::select(FullId, Id, Time, Longitude, Latitude, Depth, HPE, RMSE)

trans_IDs <- unique(vps_results$FullId)

vps_stations <- read_csv("Data/SouthAtlantic/VPS_results/stations.csv") %>% 
  select(Name = StationName, Serial, Latitude, Longitude, Depth, Start, End)


#### Load the camera data ####

camera_dat_all <- read_xlsx("Data/SouthAtlantic/SERFS_video_data_paired_sampling_2024_CTD.xlsx") %>% 
  filter(`A Video Readable` == "Yes") %>% 
  select(Project, Year, date = Date, Start_Time_GMT, 
         Station_ID, Longitude = Start_Longitude, Latitude = Start_Latitude,
         current_dir = `A Current Direction`, all_of(as.character(1:41))) %>% 
  mutate(time = as.numeric(difftime(Start_Time_GMT, as_datetime("1899-12-31 00:00:00"), units = "sec"))) %>% 
  mutate(Station_ID = paste0("Station_", row_number()))

stations <- camera_dat_all %>% 
  select(Project, Year, Station_ID, date, Longitude, Latitude, Start_Time_GMT, current_dir, time)

camera_dat <- camera_dat_all %>% 
  select(Station_ID, date, all_of(as.character(1:41))) %>% 
  pivot_longer(cols = all_of(as.character(1:41))) %>% 
  rename(frame = name, count = value) %>% 
  mutate(frame = as.numeric(frame))
camera_dat$count[is.na(camera_dat$count)] <- 0


#### Load the ROV data ####

rov_dat_raw <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/ROV/RS_2023_2024_ROV_counts_NC_FL_VPS_studies.xlsx",
                         skip = 1) %>% 
  filter(Vessel == "Jodie Lynn 2") %>% 
  mutate(Longitude = -as.numeric(conv_unit(substr(Longitude, 1, 9), from = "deg_dec_min", to = "dec_deg")),
         Latitude  = as.numeric(conv_unit(substr(Latitude, 1, 9), from = "deg_dec_min", to = "dec_deg")))

rov_dat <- bind_rows(
  rov_dat_raw %>%
    select(
      Site = `Site Number`,
      Date,
      Latitude,
      Longitude,
      area = `Transect 1 Area`,
      count = `Transect 1 RS count`
    ) %>%
    mutate(tID = 1),
  rov_dat_raw %>%
    select(
      Site = `Site Number`,
      Date,
      Latitude,
      Longitude,
      area = `Transect 2 Area`,
      count = `Transect 2 RS count`
    ) %>%
    mutate(tID = 2)
) %>%
  mutate(density = count/area)


#### ROV visualization ####

rov_summarized <- rov_dat %>%
  group_by(Longitude, Latitude) %>%
  summarize(density = mean(density))

# Map of ROV counts
rov_map <- rov_summarized %>%
  ggplot() +
  geom_point(aes(Longitude, Latitude, size = density, shape = density == 0)) +
  theme_minimal() +
  scale_color_viridis_c(end = 0.9) +
  scale_size_continuous("Snapper per m2") +
  scale_shape_manual(values = c(19, 1))

ggsave("plots/FL/rov_map.jpg", rov_map,
       width = 6, height = 5, dpi = 300)



# #### VPS visualizations ####

vps_stations %>%
  ggplot() +
  geom_point(aes(Longitude, Latitude)) +
  theme_minimal()

plot_all_fish_heatmap <- vps_results %>%
  mutate(time_abs = as.numeric(difftime(Time, min(Time), units = "hour"))) %>%
  ggplot(aes(Longitude, Latitude)) +
  geom_bin_2d() +
  geom_point(data = vps_stations, col = "darkorange") +
  theme_minimal() +
  scale_fill_viridis_c("Num. fixes", trans = "log", breaks = c(1, 20, 400, 8000), end = 0.8) +
  ggtitle("All fixes for all fish")

ggsave(plot_all_fish_heatmap, filename = "plots/FL/VPS_heatmap.jpg",
       width = 5, height = 4, dpi = 300)

plot_data_per_fish <- vps_results %>%
  count(FullId) %>%
  ggplot() +
  geom_histogram(aes(n)) +
  scale_x_continuous(trans = "log", breaks = 5^(0:6)) +
  xlab("Number of fixes per tagged fish") +
  ylab("Number of fish") +
  theme_minimal()

ggsave(plot_data_per_fish, filename = "plots/FL/VPS_perfish.jpg",
       width = 4, height = 3, dpi = 300)

vps_results %>%
  count(FullId) %>% 
  nrow()


#### Camera visualizations ####

camera_counts <- camera_dat %>%
  filter(!is.na(Station_ID)) %>%
  group_by(Station_ID, date) %>%
  summarize(total = sum(count)) %>%
  left_join(stations, by = c("date", "Station_ID")) %>%
  mutate(current_dir = ifelse(grepl("toward", tolower(current_dir)),
                              "Toward", ifelse(grepl("away", tolower(current_dir)),
                                               "Away", "Perpendicular")))

current_plot <- camera_counts %>%
  ggplot() +
  geom_point(aes(Longitude, Latitude, size = total,
                 col = current_dir, shape = total == 0)) +
  theme_bw() +
  theme(axis.ticks = element_blank()) +
  facet_wrap(~date, ncol = 1) +
  scale_color_viridis_d("Current dir.", end = 0.9) +
  scale_shape_manual("Zero count", values = c(19, 4)) +
  scale_size_continuous("Total snapper count")

time_plot <- camera_counts %>%
  ggplot() +
  geom_point(aes(Longitude, Latitude, size = total,
                 col = time / (60^2))) +
  theme_bw() +
  theme(axis.ticks = element_blank()) +
  facet_wrap(~date, ncol = 1) +
  scale_color_viridis_c(end = 0.9) +
  scale_shape_manual("Zero count", values = c(19, 4)) +
  scale_size_continuous("Total snapper count", breaks = c(1, 500, 1000))

# camera_vis <- arrangeGrob(current_plot, time_plot, nrow = 1)
ggsave("plots/FL/camera_counts.jpg", current_plot, width = 7, height = 6,
       dpi = 300)

### Map with depth and camera locations
camera_pts <- camera_counts %>%
  filter(date == as_date("2023-08-22")) %>%
  vect(, geom = c("Longitude", "Latitude"),
       crs = "+proj=longlat")

#### Camera x VPS ####

camera_summaries <- camera_counts %>%
  group_by(Longitude, Latitude) %>%
  summarize(total = sum(total))

plot_camera_x_vps <- vps_results %>%
  # filter(FullId == trans_IDs[i]) %>%
  mutate(time_abs = as.numeric(difftime(Time, min(Time), units = "hour"))) %>%
  ggplot(aes(Longitude, Latitude)) +
  geom_bin_2d() +
  # geom_point(data = vps_stations, col = "#222222", pch = 3) +
  geom_point(data = camera_summaries, aes(size = total), col = "black") +
  theme_minimal() +
  scale_size_continuous("Camera counts") +
  scale_fill_viridis_c(option = "plasma", "Num. VPS fixes",
                       trans = "log", breaks = c(1, 20, 400, 8000),
                       begin = 0.1, end = 0.9) +
  ggtitle("VPS x camera data")

ggsave("plots/FL/VPS_camera_overlay.jpg", plot_camera_x_vps,
       width = 5, height = 4, dpi = 300)



#### ROV x VPS ####

plot_camera_x_rov <- vps_results %>%
  # filter(FullId == trans_IDs[i]) %>%
  mutate(time_abs = as.numeric(difftime(Time, min(Time), units = "hour"))) %>%
  ggplot(aes(Longitude, Latitude)) +
  geom_bin_2d() +
  # geom_point(data = vps_stations, col = "#222222", pch = 3) +
  geom_point(data = rov_summarized, aes(size = density, shape = density == 0), col = "black") +
  theme_minimal() +
  scale_size_continuous("ROV density") +
  scale_fill_viridis_c(option = "plasma", "Num. VPS fixes",
                       trans = "log", breaks = c(1, 20, 400, 8000),
                       begin = 0.1, end = 0.9) +
  ggtitle("ROV x camera data") +
  scale_shape_manual(values=c(19,1))

ggsave("plots/FL/VPS_ROV_overlay.jpg", plot_camera_x_rov,
       width = 5, height = 4, dpi = 300)

