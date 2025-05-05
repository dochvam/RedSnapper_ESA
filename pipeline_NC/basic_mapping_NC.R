library(readxl)
library(tidyverse)
library(gridExtra)
library(terra)
library(tidyterra)

set.seed(5249938)

#### Load the VPS data ####

vps_results <- list.files("Data/SnapperMvmtAbundanceStudy/VPS_Data/VPS-ChickenRock-01-Results-20240202/results/animal",
                          full.names = TRUE) %>% 
  lapply(read_csv, progress = F, col_types = list(Transmitter = col_character())) %>% 
  bind_rows() %>% 
  select(FullId, Id, Time, Longitude, Latitude, Depth, HPE, RMSE)

trans_IDs <- unique(vps_results$FullId)

vps_stations <- read_xlsx("Data/SnapperMvmtAbundanceStudy/VPS_Data/2023 American Red Snapper VPS Spec.xlsx",
                          sheet = "Stations", skip = 2, n_max = 23) %>% 
  select(Name, Latitude, Longitude, Depth, Start, SO, End, EO, Height)


#### Load the camera data ####

camera_dat <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/cameratraps/RS_2023_full_reads_all_three_cameratrap_dates.xlsx") %>% 
  mutate(time = as.numeric(difftime(time, as_datetime("1899-12-31 00:00:00"), units = "sec"))) %>% 
  filter(is.na(comments) | grepl("tagged", comments)) # Filter out descent/retrieval frames

stations <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/cameratraps/RS_2023_full_reads_all_three_cameratrap_dates.xlsx", sheet = "StationData") %>% 
  mutate(Start_Time_GMT = as.numeric(difftime(Start_Time_GMT, as_datetime("1899-12-31 00:00:00"), units = "hours")))


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
  mutate(density = count/area)

#### Load the hardbottom raster ####

hb <- rast("Data/Hardbottom/ProcessedRasters/r_TNC.tif")

receiver_pts <- vect(vps_stations, geom = c("Longitude", "Latitude"),
                     crs = "+proj=longlat") %>% 
  project(crs(hb))

hb_local <- crop(hb, ext(receiver_pts))


#### Load the Chicken Rock map ####

crm <- rast("Data/Chicken_Rock_Map/NF-17-06-BH_Bachlers_Box_1m.bag")

#### ROV visualization ####

rov_summarized <- rov_dat %>% 
  group_by(Longitude, Latitude) %>% 
  summarize(density = mean(density), structured_pct = mean(structured_pct))

# Map of ROV counts
rov_map <- rov_summarized %>% 
  ggplot() +
  geom_point(aes(Longitude, Latitude, size = density,
                 col = structured_pct)) +
  theme_minimal() +
  scale_color_viridis_c(end = 0.9) + 
  scale_size_continuous("Snapper per m2")

ggsave("plots/rov_map.jpg", rov_map,
       width = 6, height = 5, dpi = 300)


# ROV counts as a fn of habitat
structured_plot <- rov_dat %>% 
  ggplot() + 
  geom_point(aes(structured_pct, density)) + 
  theme_minimal() + 
  xlab("Pct. structured habitat")

complexity_plot <- rov_dat %>% 
  ggplot() + 
  geom_point(aes(as.factor(hab_complexity), density),
             position = position_jitter(width = 0.2)) + 
  theme_minimal() + 
  xlab("Habitat complexity")


rov_habstruct_plot <- arrangeGrob(
  structured_plot, complexity_plot, nrow = 1, widths = c(1.3, 1)
)
ggsave("plots/rov_hab_struct.jpg", rov_habstruct_plot,
       width = 7, height = 3.5, dpi = 300)

#### VPS visualizations ####

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

ggsave(plot_all_fish_heatmap, filename = "plots/VPS_heatmap.jpg", 
       width = 5, height = 4, dpi = 300)

plot_data_per_fish <- vps_results %>% 
  count(FullId) %>% 
  ggplot() + 
  geom_histogram(aes(n)) +
  scale_x_continuous(trans = "log", breaks = 5^(0:6)) +
  xlab("Number of fixes per tagged fish") +
  ylab("Number of fish") +
  theme_minimal()

ggsave(plot_data_per_fish, filename = "plots/VPS_perfish.jpg", 
       width = 4, height = 3, dpi = 300)


#### Camera visualizations ####

camera_counts <- camera_dat %>% 
  filter(!is.na(Station_ID)) %>% 
  group_by(Station_ID, date) %>% 
  summarize(total = sum(total)) %>% 
  left_join(stations, by = c("date" = "Date", "Station_ID")) %>%
  rename(Longitude = Start_Longitude, Latitude = Start_Latitude) %>% 
  mutate(current_dir = ifelse(grepl("toward", tolower(`Current Direction`)),
                              "Toward", ifelse(grepl("away", tolower(`Current Direction`)),
                                               "Away", "Perpendicular"))) %>% 
  mutate(Longitude = as.numeric(Longitude))

current_plot <- camera_counts %>% 
  ggplot() + 
  geom_point(aes(Longitude, Latitude, size = total,
                 col = current_dir, shape = total == 0)) +
  theme_bw() +
  theme(axis.ticks = element_blank()) +
  facet_wrap(~date) +
  scale_color_viridis_d(end = 0.9) + 
  scale_shape_manual("Zero count", values = c(19, 4)) +
  scale_size_continuous("Total snapper count", breaks = c(1, 500, 1000))

time_plot <- camera_counts %>%
  ggplot() +
  geom_point(aes(Longitude, Latitude, size = total,
                 col = Start_Time_GMT)) +
  theme_bw() +
  theme(axis.ticks = element_blank()) +
  facet_wrap(~date) +
  scale_color_viridis_c(end = 0.9) +
  scale_shape_manual("Zero count", values = c(19, 4)) +
  scale_size_continuous("Total snapper count", breaks = c(1, 500, 1000))

# camera_vis <- arrangeGrob(current_plot, time_plot, nrow = 1)
ggsave("plots/camera_counts.jpg", current_plot, width = 11, height = 4,
       dpi = 300)

### Map with depth and camera locations
camera_pts <- camera_counts %>% 
  filter(date == as_date("2023-08-22")) %>% 
  vect(, geom = c("Longitude", "Latitude"),
       crs = "+proj=longlat")

depth_map <- ggplot() + 
  geom_spatraster(data = crm, aes(fill = Elevation)) +
  geom_spatvector(data = camera_pts, aes(size = total, shape = total == 0)) +
  theme_bw() +
  theme(axis.ticks = element_blank()) +
  scale_shape_manual("Zero count", values = c(19, 4)) +
  scale_size_continuous("Total snapper count", breaks = c(1, 500, 1000)) +
  scale_fill_viridis_c("Depth", begin = 0.2, end = 0.99)

hb_map <- ggplot() + 
  geom_spatraster(data = hb_local, show.legend = F) +
  geom_spatvector(data = camera_pts, aes(size = total, shape = total == 0)) +
  theme_bw() +
  theme(axis.ticks = element_blank()) +
  scale_shape_manual("Zero count", values = c(19, 4)) +
  scale_size_continuous("Total snapper count", breaks = c(1, 500, 1000))


#### Camera x VPS ####

camera_summaries <- camera_counts %>% 
  group_by(Longitude, Latitude) %>% 
  summarize(total = sum(total))

plot_camera_x_vps <- vps_results %>% 
  # filter(FullId == trans_IDs[i]) %>% 
  mutate(time_abs = as.numeric(difftime(Time, min(Time), units = "hour"))) %>% 
  ggplot(aes(Longitude, Latitude)) +
  geom_bin_2d() +
  geom_point(data = vps_stations, col = "#222222", pch = 3) +
  geom_point(data = camera_summaries, aes(size = total), col = "black") +
  theme_minimal() +
  scale_size_continuous("Camera counts") +
  scale_fill_viridis_c(option = "plasma", "Num. VPS fixes", 
                       trans = "log", breaks = c(1, 20, 400, 8000), 
                       begin = 0.1, end = 0.9) +
  ggtitle("VPS x camera data") 

ggsave("plots/VPS_camera_overlay.jpg", plot_camera_x_vps, 
       width = 5, height = 4, dpi = 300)



#### ROV x VPS ####

plot_camera_x_rov <- vps_results %>% 
  # filter(FullId == trans_IDs[i]) %>% 
  mutate(time_abs = as.numeric(difftime(Time, min(Time), units = "hour"))) %>% 
  ggplot(aes(Longitude, Latitude)) +
  geom_bin_2d(alpha = 0.7) +
  geom_point(data = vps_stations, col = "#222222", pch = 3) +
  geom_point(data = rov_summarized, aes(size = density), col = "black") +
  theme_minimal() +
  scale_size_continuous("RS count density on ROVs") +
  scale_fill_viridis_c(option = "plasma", "Num. VPS fixes", 
                       trans = "log", breaks = c(1, 20, 400, 8000), 
                       begin = 0.1, end = 0.9) +
  ggtitle("ROV x VPS data")

ggsave("plots/VPS_ROV_overlay.jpg", plot_camera_x_rov, 
       width = 5, height = 4, dpi = 300)


#### VPS x Hard bottom raster ####

hardbottom <- rast("Data/Chicken_Rock_Map/ChickenRock_Classification.tif")

vps_pts <- vps_results %>% 
  # sample_frac(0.1) %>% 
  vect(geom = c("Longitude", "Latitude"),
                crs = "+proj=longlat") %>% 
  project(crs(hardbottom))

hbmap <- ggplot() +
  geom_spatraster(data = hardbottom) +
  theme_minimal() +
  scale_fill_viridis_d() + ggtitle("Hardbottom map")

hbmap_wvps <- ggplot() +
  geom_spatraster(data = hardbottom) +
  geom_spatvector(data = vps_pts, alpha = 0.05, col = "gray") +
  theme_minimal() +
  scale_fill_viridis_d() + ggtitle("Hardbototm map w/ VPS fix overlay")

hb_all <- arrangeGrob(hbmap, hbmap_wvps, nrow = 2)
ggsave("plots/hardbottom_map.jpg", hb_all, width = 6, height = 6)


camera_pts <- camera_counts %>% 
  # filter(date == as_date("2023-08-22")) %>% 
  vect(, geom = c("Longitude", "Latitude"),
       crs = "+proj=longlat")

hbmap_wcams <- ggplot() +
  geom_spatraster(data = hardbottom) +
  geom_spatvector(data = camera_pts, col = "red", size = 1) +
  theme_minimal() +
  scale_fill_viridis_d() + 
  ggtitle("Hardbottom map (w camera locs)")
ggsave("plots/hardbottom_map_wcam.jpg", hbmap_wcams, width = 6, height = 3)





#### Plot the hardbottom mask w cameras ####
hardbottom <- rast("Data/Chicken_Rock_Map/ChickenRock_Classification.tif")

terra::values(hardbottom) <- ifelse(terra::values(hardbottom) == 4, 1, 0)
terra::values(hardbottom)[is.na(terra::values(hardbottom))] <- 0

hb_mask <- hardbottom %>% 
  aggregate(30, fun = "max")
terra::values(hb_mask)[is.nan(terra::values(hb_mask))] <- 0

hbmask_wcams <- ggplot() +
  geom_spatraster(data = hb_mask) +
  geom_spatvector(data = camera_pts, col = "black", size = 2) +
  theme_minimal() +
  scale_fill_viridis_c(begin = 0.1, end = 0.8) + 
  ggtitle("Hardbottom 30 m mask (w camera locs)")
ggsave("plots/hardbottom_modelmask_wcam.jpg", hbmask_wcams, width = 6, height = 3)




hb_raw <- rast("Data/Chicken_Rock_Map/ChickenRock_Classification.tif")

# Get all the VPS fixes of living fish
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
  project(crs(hb_raw))

e <- ext(fish_pts)
xmin(e) <- floor(xmin(e) / 50) * 50
xmax(e) <- ceiling(xmax(e) / 50) * 50
ymin(e) <- floor(ymin(e) / 50) * 50
ymax(e) <- ceiling(ymax(e) / 50) * 50

template_grid <- rast(e, res = 50)
terra::values(template_grid) <- 1:ncell(template_grid)


cell_counts <- count(extract(template_grid, fish_pts), lyr.1)
cell_counts <- left_join(data.frame(lyr.1 = 1:ncell(template_grid)), cell_counts)
cell_counts$n[is.na(cell_counts$n)] <- 0
cell_counts$prob <- (cell_counts$n+1) / sum(cell_counts$n+1)

vps_intensity_ras <- template_grid
terra::values(vps_intensity_ras) <- cell_counts$prob
crs(vps_intensity_ras) <- crs(fish_pts)


vpssurface_wcams <- ggplot() +
  geom_spatraster(data = vps_intensity_ras) +
  geom_spatvector(data = camera_pts, col = "black", size = 2) +
  theme_minimal() +
  scale_fill_viridis_c(begin = 0.1, end = 0.8) + 
  ggtitle("VPS intensity surface 30 m (w camera locs)")
ggsave("plots/VPS_surface.jpg", vpssurface_wcams, width = 6, height = 3)

terra::values(vps_intensity_ras) <- cell_counts$n+1

vpssurface_log <- ggplot() +
  geom_spatraster(data = vps_intensity_ras) +
  # geom_spatvector(data = camera_pts, col = "black", size = 2) +
  theme_minimal() +
  scale_fill_viridis_c(begin = 0.1, end = 0.8, trans = "log",
                       breaks = c(1, 20, 400, 8000)) + 
  ggtitle("VPS intensity surface 30 m (w camera locs), log-scale") +
  theme_void()
ggsave("plots/VPS_surface_logscale.jpg", vpssurface_log, width = 6, height = 3)

vpssurface_log <- ggplot() +
  geom_spatraster(data = vps_intensity_ras) +
  geom_spatvector(data = camera_pts, col = "black", aes(size = total)) +
  theme_minimal() +
  scale_fill_viridis_c(begin = 0.1, end = 0.8, trans = "log",
                       breaks = c(1, 20, 400, 8000)) + 
  ggtitle("VPS intensity surface 30 m (w camera locs), log-scale") +
  facet_wrap(~date) + 
  theme_void()
ggsave("plots/VPS_surface_wcam.jpg", vpssurface_log, width = 12, height = 3)


#### plot the hardbottom mask with ROV transects ####
ROV_transect_locs <- read_csv("Data/SnapperMvmtAbundanceStudy/CountData/ROV/ROV_latlong_real_2023_processed.csv")

df_split <- split(ROV_transect_locs, 
                  paste0(ROV_transect_locs$Site_ID, ROV_transect_locs$Date, ROV_transect_locs$transect))
line_geoms <- lapply(df_split, 
                     function(group) {
  as.matrix(group[, c("longitude", "latitude")])
})

# Create the line SpatVector
lines_vect <- vect(line_geoms, type = "lines", crs = "+proj=longlat") %>% 
  project(crs(vps_intensity_ras))

line_attributes <- do.call(rbind, lapply(df_split, function(group) {
  group[1, c("Site_ID", "Date", "transect")]  # retain one row per group
}))

# Assign attributes to the spatial object
values(lines_vect) <- line_attributes


vpssurface_wROVs_log <- ggplot() +
  geom_spatraster(data = vps_intensity_ras) +
  geom_spatvector(data = lines_vect, col = "black", linewidth = 0.75) +
  # geom_label(data = as.data.frame(terra::centroids(lines_vect), geom = "XY"),
  #            aes(x = x, y = y, label = Site_ID),
  #            size = 2,
  #            fontface = "bold",
  #            color = "black") +
  theme_minimal() +
  scale_fill_viridis_c(begin = 0.1, end = 0.8, trans = "log") + 
  ggtitle("VPS intensity surface 30 m (w camera locs), log-scale")+
  facet_wrap(~Date)

ggsave("plots/ROV_real_transects_overlay.jpg", vpssurface_wROVs_log, width = 6, height = 3)




# Check relationship between reported and apparent lengths
rov_dat_raw <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/ROV/ROV_NC_2023_Formatted_.xlsx")

rov_dat <- bind_rows(
  rov_dat_raw %>% 
    select(
      Site_ID = `Site:`,
      Date = `Date:`,
      Latitude = LAT, 
      Longitude = LONG,
      length = `Transect Length1`,
      area = `Area Surveyed1`,
      hab_complexity = `Habitat Complexity1`,
      structured_pct = `Structured Habitat %Cover1`,
      count = Count1
    ) %>% 
    mutate(Transect = 1),
  rov_dat_raw %>% 
    select(
      Site_ID = `Site:`,
      Date = `Date:`,
      Latitude = LAT, 
      Longitude = LONG,
      length = `Transect Length2`,
      area = `Area Surveyed2`,
      hab_complexity = `Habitat Complexity2`,
      structured_pct = `Structured Habitat %Cover2`,
      count = count2
    ) %>% 
    mutate(Transect = 2)
) %>% 
  select(Site_ID, Transect, Date, length_reported = length)

lines_df <- as.data.frame(lines_vect)
lines_df$length_ll <- perim(lines_vect)

left_join(rov_dat, lines_df, by = c("Site_ID", "Transect" = 'transect', "Date")) %>% 
  ggplot() + 
  geom_point(aes(length_reported, length_ll)) +
  geom_abline(slope = 1, intercept = 0)
