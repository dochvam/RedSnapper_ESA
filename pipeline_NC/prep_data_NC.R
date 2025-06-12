

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




vps_baseline <- 0.01

# Count the number of VPS fixes in each grid cell, then get a total proportion
cell_counts <- count(extract(template_grid, fish_pts), lyr.1)
cell_counts <- left_join(data.frame(lyr.1 = 1:ncell(template_grid)), cell_counts)
cell_counts$n[is.na(cell_counts$n)] <- 0
# Add 1 obs. to every cell so we don't make it impossible for a fish to be somewhere
cell_counts$prob <- (cell_counts$n+vps_baseline) / sum(cell_counts$n+vps_baseline)

# Formatting stuff
covariate_ras <- template_grid
terra::values(covariate_ras) <- cell_counts$prob
plot(covariate_ras)

covariate_mtx <- t(as.matrix(covariate_ras, wide = T))
covariate_mtx <- covariate_mtx[, ncol(covariate_mtx):1]

grid_bbox <- ext(covariate_ras)
grid_bbox <- as.numeric(c(grid_bbox[1], grid_bbox[2], grid_bbox[3], grid_bbox[4]))

x_offset <- mean(grid_bbox[1:2])
y_offset <- mean(grid_bbox[3:4])

grid_bbox[1:2] <- grid_bbox[1:2] - x_offset
grid_bbox[3:4] <- grid_bbox[3:4] - y_offset

resoln <- res(covariate_ras)[1]

# Get the informed prior for sigma
log_sigma_estimate <- read_csv("pipeline_NC/NC_results/log_sigma_estimate_NC.csv")

if (sigma_type == "Prior_Variability") {
  log_sigma_estimate <- log_sigma_estimate %>% 
    filter(type == "variability", t_interval == max(t_interval))
} else {
  log_sigma_estimate <- log_sigma_estimate %>% 
    filter(type == "mean", t_interval == max(t_interval))
}

attraction_dist_sigma <- 27.11

# log_sigma_point_estimate <- log(sqrt(exp(3.49)^2 + attraction_dist_sigma^2))

# if (sigma_type == "Mean") {
#   log_sigma_estimate <- log_sigma_estimate %>% filter(type == "mean")
# } else {
#   log_sigma_estimate <- log_sigma_estimate %>% filter(type == "variability")
# }

#### Load both days of ROV data ####

surface_to_buffer <- covariate_ras
source("uSCR_real/process_ROV_buffers.R")

processed_buffers_result <- process_ROV_buffers(buffer_size = 4.9, surface_to_buffer = surface_to_buffer)
rbe <- processed_buffers_result$rbe
rbs <- processed_buffers_result$rbs
intersections_df <- processed_buffers_result$intersections_df
rov_dat <- processed_buffers_result$rov_dat

#### Load and format the camera data ####

# Get the real camera locations and deployment times

suppressWarnings(
  camera_dat <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/cameratraps/RS_2023_full_reads_all_three_cameratrap_dates.xlsx") %>% 
    mutate(time = as.numeric(difftime(time, as_datetime("1899-12-31 00:00:00"), units = "sec"))) %>% 
    filter(camera == "A", time > 60 * 10, time <= 60 * 30) # Filter out descent/retrieval frames
)

camera_detmtx <- camera_dat %>%
  dplyr::select(Station_ID, date, total) %>% 
  group_by(Station_ID, date) %>% 
  mutate(col = paste0("V", row_number())) %>% 
  pivot_wider(names_from = col, values_from = total)
obs_cols <- colnames(camera_detmtx)[3:ncol(camera_detmtx)]
camera_detmtx$date_index <- as.numeric(as.factor(camera_detmtx$date))

if (exists("data_thin_interval")) {
  obs_cols <- obs_cols[1 + 0:(length(obs_cols)/data_thin_interval-1) * data_thin_interval]
}

camera_locations_corrected <- read_csv("intermediate/corrected_camera_stations.csv") %>% 
  dplyr::select(Station_ID, Date, Longitude, Latitude)

suppressWarnings(
  camera_locs <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/cameratraps/RS_2023_full_reads_all_three_cameratrap_dates.xlsx", sheet = "StationData") %>% 
    filter(`Camera (A or C)` == "A") %>%
    dplyr::select(Station_ID, Date, Time = Start_Time_GMT, 
                  current_dir = `Current Direction`) %>% 
    left_join(camera_locations_corrected, by = c("Station_ID", "Date"))
)

camera_locs$current_dir <- ifelse(
  grepl("towards", tolower(camera_locs$current_dir)), 1, 
  ifelse(grepl("away", tolower(camera_locs$current_dir)), 2, 3)
)
camera_locs$date_index <- as.numeric(as.factor(camera_locs$Date))

if (exists("target_date")) {
  camera_detmtx <- camera_detmtx %>% filter(date_index == target_date)
  camera_locs <- camera_locs %>% filter(date_index == target_date)
}


y_mtx <- as.matrix(camera_detmtx[, obs_cols])

# Make sure order of camera locs matches order of detection history
stopifnot(all(camera_locs$Station_ID == camera_detmtx$Station_ID) &
            all(camera_locs$Date == camera_detmtx$date))

camera_locs$Longitude <- as.numeric(camera_locs$Longitude)

camera_pts <- vect(camera_locs, geom = c("Longitude", "Latitude"),
                   crs = "+proj=longlat") %>% 
  project(crs(hb_raw)) %>% 
  as.data.frame(geom = "XY")

X <- camera_pts[, c("x", "y")]
X[, 1] <- X[, 1] - x_offset
X[, 2] <- X[, 2] - y_offset

# For each camera, get "time since first deployment began"
camera_locs <- camera_locs %>% 
  group_by(Date) %>% 
  mutate(mins_since_start = as.numeric(difftime(Time, min(Time), units = "mins")) + 20)

Dates <- unique(camera_locs$Date)

#### Data reformatting for model ####


if (!exists("target_date")) {
  # Make a 3D array of camera coords  
  mtx_nrow <- max(table(camera_locs$Date))
  ncam_vec <- numeric(3)
  X_array <- array(dim = c(mtx_nrow, 3, 2))
  for (i in 1:3) {
    ncam_vec[i] <- sum(camera_locs$Date == Dates[i])
    X_array[1:ncam_vec[i], i, 1] <- X$x[camera_locs$Date == Dates[i]]
    X_array[1:ncam_vec[i], i, 2] <- X$y[camera_locs$Date == Dates[i]]
  }
  # Convert to a 2d matrix
  X_mtx <- rbind(X_array[,1,], X_array[,2,], X_array[1:14,3,])
  X_datevec <- c(rep(1, 18), rep(2, 18), rep(3, 14))
  stopifnot(nrow(X_mtx) == length(X_datevec))
  
  idate <- rep(1:3, each = M/3)
} else {
  X_mtx <- X
  X_datevec <- rep(1, nrow(X_mtx))
  
  idate <- rep(1, M)
}

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
                  spatial_type = spatial_type,
                  log_sigma_prior_mean = log_sigma_estimate$mean,
                  log_sigma_prior_sd = log_sigma_estimate$sd,
                  current_dir = camera_locs$current_dir,
                  K1D = rep(K, J),
                  n.samples = n.samples,
                  xlim = xlim,
                  ylim = ylim,
                  datevec = X_datevec,
                  idate = idate,
                  hm = covariate_mtx,
                  ones_mtx = ones_mtx,
                  hm_nrow = nrow(covariate_mtx),
                  hm_ncol = ncol(covariate_mtx),
                  resoln = res(covariate_ras)[1],
                  nROV = nrow(rov_dat),
                  rb_weights = intersections_df$weight,
                  rov_cell_xvec = intersections_df$x_ind,
                  rov_cell_yvec = intersections_df$y_ind,
                  rbe = rbe, rbs = rbs,
                  sigma_type = sigma_type,
                  attraction_dist_sigma = attraction_dist_sigma
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
log_sigma_init <- constants$log_sigma_prior_mean
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
                 ROV_offset = 1,
                 spatial_beta = 0.3)

if (sigma_type != "Constant") {
  Niminits[["log_sigma_long"]] <- log_sigma_init
}

Nimdata <- list(y.true=matrix(NA,nrow=(M),ncol=J),
                ROV_obs = rov_dat$count,
                ID = rep(NA, n.samples),
                z = rep(NA, M),
                X = as.matrix(X_mtx),
                capcounts=rep(NA, M))
