library(tidyverse)
library(terra)
library(readxl)
library(measurements)

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


vps_baseline <- 0.01

if (spatial_type == "VPS_map") {
  # Count the number of VPS fixes in each grid cell, then get a total proportion
  cell_counts <- count(extract(template_grid, fish_pts), lyr.1)
  cell_counts <- left_join(data.frame(lyr.1 = 1:ncell(template_grid)), cell_counts)
  cell_counts$n[is.na(cell_counts$n)] <- 0
  # Add 1 obs. to every cell so we don't make it impossible for a fish to be somewhere
  cell_counts$prob <- (cell_counts$n+vps_baseline) / sum(cell_counts$n+vps_baseline)
  
  # Formatting stuff
  covariate_ras <- template_grid
  terra::values(covariate_ras) <- cell_counts$prob + vps_baseline

  covariate_mtx <- t(as.matrix(covariate_ras, wide = T))
  covariate_mtx <- covariate_mtx[, ncol(covariate_mtx):1]
  
} else if (spatial_type == "HB_map") {
  hbmap <- vect("Data/SouthAtlantic/GraphicsLayer_G2F1.shp") %>% 
    project("ESRI:102003")
  
  covariate_ras <- rasterize(hbmap, template_grid, touches = T)
  terra::values(covariate_ras)[is.nan(terra::values(covariate_ras))] <- vps_baseline
  
  covariate_mtx <- t(as.matrix(covariate_ras, wide = T))
  covariate_mtx <- covariate_mtx[, ncol(covariate_mtx):1]
}

grid_bbox <- ext(covariate_ras)
grid_bbox <- as.numeric(c(grid_bbox[1], grid_bbox[2], grid_bbox[3], grid_bbox[4]))

x_offset <- mean(grid_bbox[1:2])
y_offset <- mean(grid_bbox[3:4])

grid_bbox[1:2] <- grid_bbox[1:2] - x_offset
grid_bbox[3:4] <- grid_bbox[3:4] - y_offset

resoln <- res(covariate_ras)[1]

# Get the informed prior for sigma
log_sigma_estimate <- read_csv("pipeline_FL/FL_results/log_sigma_estimate_FL.csv")
if (sigma_type == "Mean") {
  log_sigma_estimate <- log_sigma_estimate %>% filter(type == "mean")
} else {
  log_sigma_estimate <- log_sigma_estimate %>% filter(type == "variability")
}


#### Load both days of ROV data ####

surface_to_buffer <- covariate_ras

# process ROV buffers. The *avg* density is determined by the covariate levels
#       in a 100 m buffer around the ROV point, but the area sampled is determined
#       by the reported area, which can be calculated as a fraction of the total
#       sampling area (so E = pct. of survey area * pct of fish * N/2)

rov_dat_raw <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/ROV/RS_2023_2024_ROV_counts_NC_FL_VPS_studies.xlsx",
                         skip = 1) %>% 
  filter(Vessel == "Jodie Lynn 2") %>% 
  mutate(Longitude = -as.numeric(conv_unit(substr(Longitude, 1, 9), from = "deg_dec_min", to = "dec_deg")),
         Latitude  = as.numeric(conv_unit(substr(Latitude, 1, 9), from = "deg_dec_min", to = "dec_deg")))

rov_dat <- bind_rows(
  rov_dat_raw %>%
    dplyr::select(
      Site = `Site Number`,
      Date,
      Latitude,
      Longitude,
      area = `Transect 1 Area`,
      count = `Transect 1 RS count`
    ) %>%
    mutate(tID = 1),
  rov_dat_raw %>%
    dplyr::select(
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


rov_pts <- vect(rov_dat, geom = c("Longitude", "Latitude"),
                crs = "+proj=longlat") %>% 
  project("ESRI:102003")

buffer_size <- 100

rov_dat$buffer_area <- NA
rov_dat$surface_mean <- NA
rov_dat$surface_sum <- NA

intersections_df <- data.frame()

for (i in 1:nrow(rov_dat)) {
  this_pt <- rov_pts[i, ]
  
  buffer <- buffer(this_pt, buffer_size)
  
  this_value <- terra::extract(surface_to_buffer, buffer, weights = TRUE, cells = TRUE)
  
  rov_dat$buffer_area[i]  <- expanse(buffer)
  rov_dat$surface_mean[i] <- sum(this_value$lyr.1 * (this_value$weight / sum(this_value$weight)))
  rov_dat$surface_sum[i]  <- sum(this_value$lyr.1 * this_value$weight)
  
  
  # Down-weight all weights based on the reported area covered by the transect
  area_correction <- rov_dat$area[i] / rov_dat$buffer_area[i]
  this_value$weight <- this_value$weight * area_correction
  
  this_value <- bind_cols(this_value, rov_dat[i, c("Site", "Date", "tID")])
  intersections_df <- bind_rows(intersections_df, this_value)
}

intersections_df <- intersections_df %>% 
  mutate(
    y_ind = 1 + dim(surface_to_buffer)[1] - ceiling(cell / dim(surface_to_buffer)[2]),
    x_ind = (cell %% dim(surface_to_buffer)[2]),
    ROV_ID = paste(Site, Date, tID, sep = "_")
  ) %>% 
  arrange(ROV_ID)


rov_dat <- rov_dat %>% 
  mutate(ROV_ID = paste(Site, Date, tID, sep = "_")) %>% 
  arrange(ROV_ID) %>% 
  filter(!is.na(buffer_area))

rbs <- rbe <- numeric(nrow(rov_dat))
for (i in 1:length(rbs)) {
  rbs[i] <- min(which(intersections_df$ROV_ID == rov_dat$ROV_ID[i]))
  rbe[i] <- max(which(intersections_df$ROV_ID == rov_dat$ROV_ID[i]))
}


#### Load and format the camera data ####

# Get the real camera locations and deployment times

camera_dat_all <- read_xlsx("Data/SouthAtlantic/SERFS_video_data_paired_sampling_2024_CTD.xlsx") %>% 
  filter(`A Video Readable` == "Yes") %>% 
  dplyr::select(Project, Year, date = Date, Start_Time_GMT, 
         Station_ID, Longitude = Start_Longitude, Latitude = Start_Latitude,
         current_dir = `A Current Direction`, all_of(as.character(1:41))) %>% 
  mutate(time = as.numeric(difftime(Start_Time_GMT, as_datetime("1899-12-31 00:00:00"), units = "sec"))) %>% 
  mutate(Station_ID = paste0("Station_", row_number()))

camera_locs <- camera_dat_all %>% 
  dplyr::select(Project, Year, Station_ID, date, Longitude, Latitude, Start_Time_GMT, current_dir, time)

camera_dat <- camera_dat_all %>% 
  dplyr::select(Station_ID, date, all_of(as.character(1:41))) %>% 
  pivot_longer(cols = all_of(as.character(1:41))) %>% 
  rename(frame = name, count = value) %>% 
  mutate(frame = as.numeric(frame))
camera_dat$count[is.na(camera_dat$count)] <- 0

camera_detmtx <- camera_dat %>%
  dplyr::select(Station_ID, date, total = count) %>% 
  group_by(Station_ID, date) %>% 
  mutate(col = paste0("V", row_number())) %>% 
  pivot_wider(names_from = col, values_from = total)
obs_cols <- colnames(camera_detmtx)[3:ncol(camera_detmtx)]

y_mtx <- as.matrix(camera_detmtx[, obs_cols])



camera_locs$current_dir <- ifelse(
  grepl("towards", tolower(camera_locs$current_dir)), 1, 
  ifelse(grepl("away", tolower(camera_locs$current_dir)), 2, 3)
)

# Make sure order of camera locs matches order of detection history
stopifnot(all(camera_locs$Station_ID == camera_detmtx$Station_ID) &
            all(camera_locs$date == camera_detmtx$date))

camera_locs$Longitude <- as.numeric(camera_locs$Longitude)

camera_pts <- vect(camera_locs, geom = c("Longitude", "Latitude"),
                   crs = "+proj=longlat") %>% 
  project("ESRI:102003") %>% 
  as.data.frame(geom = "XY")

X <- camera_pts[, c("x", "y")]
X[, 1] <- X[, 1] - x_offset
X[, 2] <- X[, 2] - y_offset

# For each camera, get "time since first deployment began"
camera_locs <- camera_locs %>% 
  group_by(date) %>% 
  mutate(mins_since_start = as.numeric(difftime(time, min(time), units = "mins")) + 20)

Dates <- unique(camera_locs$date)

#### Data reformatting for model ####


# Make a 3D array of camera coords  
mtx_nrow <- max(table(camera_locs$date))
ncam_vec <- numeric(2)
X_array <- array(dim = c(mtx_nrow, 2, 2))
for (i in 1:2) {
  ncam_vec[i] <- sum(camera_locs$date == Dates[i])
  X_array[1:ncam_vec[i], i, 1] <- X$x[camera_locs$date == Dates[i]]
  X_array[1:ncam_vec[i], i, 2] <- X$y[camera_locs$date == Dates[i]]
}
# Convert to a 2d matrix
X_mtx <- rbind(X_array[1:ncam_vec[1],1,], X_array[1:ncam_vec[2],2,])
X_datevec <- c(rep(1, ncam_vec[1]), rep(2, ncam_vec[2]))
stopifnot(nrow(X_mtx) == length(X_datevec))

#### Get observed counts ####

xlim <- grid_bbox[1:2]
ylim <- grid_bbox[3:4]

J <- nrow(X_mtx)
K <- ncol(y_mtx)

n.samples <- sum(y_mtx)

this.j <- this.k <- rep(NA,n.samples)
idx <- 1
for(j in 1:J){ #then traps
  for(k in 1:K){ #then occasions
    if(y_mtx[j,k]>0){ #is there at least one sample here?
      for(l in 1:y_mtx[j,k]){ #then samples
        this.j[idx] <- j
        this.k[idx] <- k
        idx <- idx+1
      }
    }
  }
}

# sim.out <- sim_counts_wHabCovar(
#   M = M, psi = true_psi, p0 = true_p0, sigma = true_sigma, 
#   X_mtx = X_mtx, X_datevec = X_datevec,
#   J = J, xlim = xlim, ylim = ylim, K = K,
#   resoln = resoln,
#   habitatMask = covariate_mtx,
#   spatial_beta = true_spatial_beta
# )

# this.j and this.k describe the samples and occasions on which each detection occurred

# Simulate ROV data #
# ROV_counts <- sim_ROV(nROV = nrow(rov_dat),
#                       rb_weights = intersections_df$weight,
#                       rov_cell_xvec = intersections_df$x_ind,
#                       rov_cell_yvec = intersections_df$y_ind,
#                       rbe = rbe, rbs = rbs, 
#                       spatial_beta = true_spatial_beta, 
#                       habitatMask = covariate_mtx,
#                       sim.out = sim.out)

ones_mtx <- covariate_mtx
ones_mtx[] <- 1

#### Make NIMBLE model input lists ####

constants <- list(M = M,
                  J = J,
                  integration_type = integration_type,
                  log_sigma_prior_mean = log_sigma_estimate$mean,
                  log_sigma_prior_sd = log_sigma_estimate$sd,
                  current_dir = camera_locs$current_dir,
                  K1D = rep(K, J),
                  n.samples = n.samples,
                  xlim = xlim,
                  ylim = ylim,
                  datevec = X_datevec,
                  idate = rep(1:2, each = M/2),
                  hm = covariate_mtx,
                  ones_mtx = ones_mtx,
                  hm_nrow = nrow(covariate_mtx),
                  hm_ncol = ncol(covariate_mtx),
                  resoln = res(covariate_ras)[1],
                  nROV = nrow(rov_dat),
                  rb_weights = intersections_df$weight,
                  rov_cell_xvec = intersections_df$x_ind,
                  rov_cell_yvec = intersections_df$y_ind,
                  rbe = rbe, rbs = rbs
                  )


z_init <- rbinom(M, 1, 0.5)
s_init <- initialize_s(
  M = M,
  xlim = xlim,
  ylim = ylim,
  resoln = resoln,
  habitatMask = covariate_mtx,
  spatial_beta = 0.5
)
# log_sigma_init <- constants$log_sigma_prior_mean
log_sigma_init <- 3 # <- need this to be high enough that there are no 0-prob cameras from random initial locs
sigma_init <- exp(log_sigma_init)
p0_init <- rep(0.8, 3)
y.true.init <- initialize_ytrue(M,
                                z_init, 
                                s_init, 
                                this.j = this.j,
                                this.k = this.k,
                                X_mtx = X_mtx,
                                current_dir = camera_locs$current_dir,
                                X_datevec = X_datevec, 
                                idate = constants$idate,
                                n.samples = n.samples,
                                sigma_init = sigma_init, p0_init = p0_init, K = K)

Niminits <- list(z = z_init,
                 N = sum(z_init), #must initialize N to be the sum of z init
                 lambda.N=sum(z_init), #initializing lambda.N to be consistent with N.init
                 s = s_init,
                 ID = y.true.init$ID,
                 capcounts = rowSums(y.true.init$ytrue2D),
                 y.true = y.true.init$ytrue2D,
                 p0 = p0_init, 
                 sigma = sigma_init,
                 log_sigma = log_sigma_init,
                 ROV_offset = 1,
                 spatial_beta = 0.3)


Nimdata <- list(y.true=matrix(NA,nrow=(M),ncol=J),
                ROV_obs = rov_dat$count,
                ID = rep(NA, n.samples),
                z = rep(NA, M),
                X = as.matrix(X_mtx),
                capcounts=rep(NA, M))


