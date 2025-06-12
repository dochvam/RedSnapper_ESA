library(tidyverse)
library(terra)
library(tidyterra)
library(geodata)
library(gridExtra)

get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#### Data figure panel for CR #### 

camera_df <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/cameratraps/RS_2023_full_reads_all_three_cameratrap_dates.xlsx", sheet = "StationData") %>% 
  filter(`Camera (A or C)` == "A") %>% 
  mutate(Start_Time_GMT = as.numeric(difftime(Start_Time_GMT, as_datetime("1899-12-31 00:00:00"), units = "hours")),
         Longitude = as.numeric(Start_Longitude), Latitude = Start_Latitude) %>% 
  select(Longitude, Latitude, Date) %>% 
  mutate(Name = paste0("Camera (", Date, ")"))

ROV_df <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/ROV/RS_2023_2024_ROV_counts_NC_FL_VPS_studies.xlsx", skip = 1) %>% 
  filter(Region == 1)%>% 
  mutate(Latitude = as.numeric(conv_unit(substr(Latitude, 1, 9), from = "deg_dec_min", to = "dec_deg"))) %>% 
  mutate(Longitude = -as.numeric(conv_unit(substr(Longitude, 1, 9), from = "deg_dec_min", to = "dec_deg"))) %>% 
    select(
      Latitude, 
      Longitude,
      Date
    ) %>% 
  mutate(Name = paste0("ROV (", Date, ")"))


pt_locs_CR <- bind_rows(camera_df, ROV_df) %>% 
  vect(geom = c("Longitude", "Latitude"),
       crs = "+proj=longlat")





#### Process the real data ####
# First, we process all the real data. This allows us to define simulation data
# dimensions, use the real VPS data (as spatial covar.) and camera locations,
# etc. We ultimately only simulate the observed counts at each camera and ROV

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

# Convert to pts
fish_pts <- vect(fish_positions, geom = c("Longitude", "Latitude"),
                 crs = "+proj=longlat") %>% 
  project(crs(hb_raw))

# Define a spatial extent encompassing all VPS locs, snapped to a 50 m grid
grid_buffer <- 300
grid_resoln <- 50

e <- ext(fish_pts)
xmin(e) <- floor((xmin(e) - grid_buffer) / grid_resoln) * grid_resoln
xmax(e) <- ceiling((xmax(e) + grid_buffer) / grid_resoln) * grid_resoln
ymin(e) <- floor((ymin(e) - grid_buffer) / grid_resoln) * grid_resoln
ymax(e) <- ceiling((ymax(e) + grid_buffer) / grid_resoln) * grid_resoln

template_grid <- rast(e, res = grid_resoln)
crs(template_grid) <- crs(fish_pts)
terra::values(template_grid) <- 1:ncell(template_grid)




vps_baseline <- 0.1

# Count the number of VPS fixes in each grid cell, then get a total proportion
cell_counts <- count(extract(template_grid, fish_pts), lyr.1)
cell_counts <- left_join(data.frame(lyr.1 = 1:ncell(template_grid)), cell_counts)
cell_counts$n[is.na(cell_counts$n)] <- 0
# Add 1 obs. to every cell so we don't make it impossible for a fish to be somewhere
cell_counts$prob <- (cell_counts$n+vps_baseline)

# Formatting stuff
covariate_ras_CR <- template_grid
terra::values(covariate_ras_CR) <- cell_counts$n

CR_o1 <- ggplot() +
  geom_spatraster(data = covariate_ras_CR) +
  geom_spatvector(data = pt_locs_CR, col = "black", size = 3) +
  # geom_spatvector(data = rov_pts, aes(shape = as.character(date)), col = "black", size = 2) +
  theme_minimal() +
  scale_fill_viridis_c("Num. fixes", begin = 0.1, end = 0.8, trans = "log",
                       breaks = c(1, 20, 400, 8000), limits = c(0.1, 22000)) + 
  ggtitle("(B) Chicken Rock") +
  theme_void() + 
  facet_wrap(~Name, nrow = 1)

names(covariate_ras_CR) <- "NumFixes"
CR_o2 <- ggplot() +
  geom_spatraster(data = covariate_ras_CR, show.legend = TRUE, aes(fill = NumFixes)) +
  geom_spatvector(data = pt_locs_CR[pt_locs_CR$Name == "Camera (2023-08-22)", ], col = "black", size = 2) +
  # geom_spatvector(data = rov_pts, aes(shape = as.character(date)), col = "black", size = 2) +
  scale_fill_viridis_c("Num. fixes", begin = 0.1, end = 0.8, trans = "log", na.value = NA,
                       breaks = c(1, 20, 400, 8000), limits = c(0.1, 22000)) + 
  ggtitle("(B) Chicken Rock") +
  theme_minimal() + 
  theme(legend.position = "right")



#### Data figure panel for TM #### 

camera_df <- read_xlsx("Data/SouthAtlantic/SERFS_video_data_paired_sampling_2024_CTD.xlsx") %>% 
  filter(`A Video Readable` == "Yes") %>% 
  dplyr::select(Project, Year, date = Date, Start_Time_GMT, 
                Station_ID, Longitude = Start_Longitude, Latitude = Start_Latitude,
                current_dir = `A Current Direction`, all_of(as.character(1:41))) %>% 
  mutate(time = as.numeric(difftime(Start_Time_GMT, as_datetime("1899-12-31 00:00:00"), units = "sec"))) %>% 
  mutate(Station_ID = paste0("Station_", row_number())) %>% 
  dplyr::select(Date = date, Longitude, Latitude) %>% 
  mutate(Name = paste0("Camera (", Date, ")"))

ROV_df <- rov_dat_raw <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/ROV/RS_2023_2024_ROV_counts_NC_FL_VPS_studies.xlsx",
                                   skip = 1) %>% 
  filter(Vessel == "Jodie Lynn 2") %>% 
  mutate(Longitude = -as.numeric(conv_unit(substr(Longitude, 1, 9), from = "deg_dec_min", to = "dec_deg")),
         Latitude  = as.numeric(conv_unit(substr(Latitude, 1, 9), from = "deg_dec_min", to = "dec_deg"))) %>% 
  select(
    Latitude, 
    Longitude,
    Date
  ) %>% 
  mutate(Name = paste0("ROV (", Date, ")"))


pt_locs_TM <- bind_rows(camera_df, ROV_df) %>% 
  vect(geom = c("Longitude", "Latitude"),
       crs = "+proj=longlat")



#### Process the real data ####
# First, we process all the real data. This allows us to define simulation data
# dimensions, use the real VPS data (as spatial covar.) and camera locations,
# etc. We ultimately only simulate the observed counts at each camera and ROV

# Get all the VPS fixes of living fish
fish_positions <- list.files("Data/SouthAtlantic/VPS_results/animal/",
                             full.names = TRUE) %>% 
  lapply(read_csv, progress = F, col_types = list(Transmitter = col_character())) %>% 
  bind_rows() %>% 
  dplyr::select(FullId, Id, Time, Longitude, Latitude, Depth, HPE, RMSE)


fish_pts <- vect(fish_positions, geom = c("Longitude", "Latitude"),
                 crs = "+proj=longlat") %>% 
  project("ESRI:102003")



# Define a spatial extent encompassing all VPS locs, snapped to a 100 m grid
grid_buffer <- 300
grid_resoln <- 50

e <- ext(fish_pts)
xmin(e) <- floor((xmin(e) - grid_buffer) / grid_resoln) * grid_resoln
xmax(e) <- ceiling((xmax(e) + grid_buffer) / grid_resoln) * grid_resoln
ymin(e) <- floor((ymin(e) - grid_buffer) / grid_resoln) * grid_resoln
ymax(e) <- ceiling((ymax(e) + grid_buffer) / grid_resoln) * grid_resoln

template_grid <- rast(e, res = grid_resoln)
crs(template_grid) <- crs(fish_pts)
terra::values(template_grid) <- 1:ncell(template_grid)


vps_baseline <- 0.1

# Count the number of VPS fixes in each grid cell, then get a total proportion
cell_counts <- count(extract(template_grid, fish_pts), lyr.1)
cell_counts <- left_join(data.frame(lyr.1 = 1:ncell(template_grid)), cell_counts)
# cell_counts$n[is.na(cell_counts$n)] <- 0
# Add 1 obs. to every cell so we don't make it impossible for a fish to be somewhere
# cell_counts$prob <- (cell_counts$n+vps_baseline)

# Formatting stuff
covariate_ras_TM <- template_grid
terra::values(covariate_ras_TM) <- cell_counts$n



TM_o1 <- ggplot() +
  geom_spatraster(data = covariate_ras_TM) +
  geom_spatvector(data = pt_locs_TM, col = "black", size = 3) +
  # geom_spatvector(data = rov_pts, aes(shape = as.character(date)), col = "black", size = 2) +
  theme_minimal() +
  scale_fill_viridis_c("Num. fixes", begin = 0.1, end = 0.8, trans = "log", na.value = "#eeeeee",
                       breaks = c(1, 20, 400, 8000), limits = c(0.1, 22000)) + 
  ggtitle("(C) Turtle Mount") +
  theme_void() + 
  facet_wrap(~Name, nrow = 1)


TM_o2 <- ggplot() +
  geom_spatraster(data = covariate_ras_TM) +
  geom_spatvector(data = pt_locs_TM[grepl("Camera", pt_locs_TM$Name), ], 
                  col = "black", size = 2) +
  # geom_spatvector(data = rov_pts, aes(shape = as.character(date)), col = "black", size = 2) +
  theme_minimal() +
  scale_fill_viridis_c("Num. fixes", begin = 0.1, end = 0.8, trans = "log", na.value = NA,#na.value = "#eeeeee",
                       breaks = c(1, 20, 400, 8000), limits = c(0.1, 22000)) + 
  ggtitle("(C) Turtle Mount") +
  theme_minimal() 

legend <- get_legend(TM_o2)


#### Inset showing location of each site on the east coast #### 
map <- gadm("USA", path = "data", level = 1, resolution = 2)

CR_location <- vect(
  data.frame(name = 'Chicken Rock', longitude = -76.14194, latitude = 34.62383),
  geom = c("longitude", "latitude"),
  crs = "+proj=longlat"
)

TM_location <- vect(
  data.frame(name = 'Turtle Mount', longitude = -80.47271, latitude = 29.02822),
  geom = c("longitude", "latitude"),
  crs = "+proj=longlat"
)

map_inset <- ggplot() +
  geom_spatvector(data = map, fill = "white") +
  geom_spatvector(data = CR_location, color = "black", size = 3) +
  geom_spatvector_label(data = CR_location, aes(label = name), 
                        color = "black", size = 2, nudge_y = 0.6) +
  geom_spatvector(data = TM_location, color = "black", size = 3) +
  geom_spatvector_label(data = TM_location, aes(label = name), 
                        color = "black", size = 2, nudge_y = 0.6) +
  theme_minimal() +  xlab("") + ylab("") +
  # theme(panel.background = element_rect(fill = "#cee8f2", linewidth = 0)) +
  ylim(25, 38) + 
  xlim(-83, -75) +
  ggtitle("(A) Map\n")


#### Put it all together ####

layout_mtx <- matrix(c(1, 1, 2, 3, 2, NA), nrow = 2)

spat_fig_o1 <- arrangeGrob(map_inset + theme(plot.margin = margin(t = 20,  # Top margin
                                                                  r = 20,  # Right margin
                                                                  b = 20,  # Bottom margin
                                                                  l = 20)), 
                           CR_o1 + theme(legend.position = "None"), 
                           TM_o1 + theme(legend.position = "None"), 
                           layout_matrix = layout_mtx,
                           widths = c(2, 5, 1.2), 
                           right = legend)
ggsave(plot = spat_fig_o1, filename = "plots/Figure2_Data_v1.jpg",
       width = 11, height = 4, dpi = 300)


# layout_mtx <- matrix(c(1, 1, 2, 3), nrow = 2)

spat_fig_o2 <- arrangeGrob(map_inset,
                           CR_o2 + theme(legend.position = "None"), 
                           TM_o2 + theme(legend.position = "None"), 
                           # layout_matrix = layout_mtx,
                           nrow = 1,
                           widths = c(2.5, 5, 4),
                           # heights = c(0.8, 1),
                           right = legend)
ggsave(plot = spat_fig_o2, filename = "plots/Figure2_Data_v2.jpg",
       width = 12, height = 5, dpi = 300)

