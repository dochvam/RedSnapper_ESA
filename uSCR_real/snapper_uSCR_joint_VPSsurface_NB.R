library(nimbleEcology)
library(tidyverse)
library(readxl)
library(terra)

#### nimbleFunction GetDetectionRate ####
GetDetectionRate <- nimbleFunction(
  run = function(s = double(1), lam0=double(1), 
                 sigma=double(0), 
                 current_dir = double(1),
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      ans <- lam0[current_dir]*exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)


dHabDistr <- nimbleFunction(run = function(x = double(1),
                                           xmax = double(0),
                                           xmin = double(0),
                                           ymax = double(0),
                                           ymin = double(0),
                                           resoln = double(0),
                                           habitatMask = double(2),
                                           log = logical(0)) {
  if (x[1] < xmin) return(-Inf)
  if (x[1] > xmax) return(-Inf)
  if (x[2] < ymin) return(-Inf)
  if (x[2] > ymax) return(-Inf)
  
  # Convert s to row/col
  row <- trunc((x[1] - xmin) / resoln) + 1
  col <- trunc((x[2] - ymin) / resoln) + 1
  
  if (log) return(log(habitatMask[row, col]))
  else return(habitatMask[row, col])
  returnType(double(0))
}, buildDerivs = FALSE)

rHabDistr <- nimbleFunction(run = function(n = double(0),
                                           xmax = double(0),
                                           xmin = double(0),
                                           ymax = double(0),
                                           ymin = double(0),
                                           resoln = double(0),
                                           habitatMask = double(2)
){
  returnType(double(1))
  
  # Randomly select a row based on their relative weights
  rowProbs <- numeric(dim(habitatMask)[1])
  for (i in 1:(dim(habitatMask)[1])) {
    rowProbs[i] <- sum(habitatMask[i, ])
  }
  row <- rcat(1, prob = rowProbs)
  
  # Randomly select a column from that row
  col <- rcat(1, prob = habitatMask[row, ])
  
  # # Generate a random coord. uniformly within that cell
  x <- c(
    xmin + resoln * (row - 1) + runif(1, 0, resoln),
    ymin + resoln * (col - 1) + runif(1, 0, resoln)
  )
  return(x)
  
}, buildDerivs = FALSE)
# ctest <- compileNimble(rHabDistr)

initialize_s <- function(n, t,
                         xmax,
                         xmin,
                         ymax,
                         ymin,
                         resoln,
                         habitatMask) {
  s <- array(NA, dim =c(n, 3, 2))
  for (i in 1:n) {
    for (ti in 1:t) {
      s[i, ti, ] <- rHabDistr(1, xmax = xmax, xmin = xmin,
                          ymax = ymax, ymin = ymin, resoln = resoln,
                          habitatMask = habitatMask)
      
    }
  }
  s
}


#### Construct the VPS intensity surface ####

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
plot(vps_intensity_ras)

vps_mtx <- t(as.matrix(vps_intensity_ras, wide = T))
vps_mtx <- vps_mtx[, ncol(vps_mtx):1]

grid_bbox <- ext(vps_intensity_ras)
grid_bbox <- as.numeric(c(grid_bbox[1], grid_bbox[2], grid_bbox[3], grid_bbox[4]))

x_offset <- mean(grid_bbox[1:2])
y_offset <- mean(grid_bbox[3:4])

grid_bbox[1:2] <- grid_bbox[1:2] - x_offset
grid_bbox[3:4] <- grid_bbox[3:4] - y_offset


#### Load the ROV data ####

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
         structured_area = area * structured_pct) %>% 
  filter(structured_area > 0)

#### Format the camera data ####
M <- 5000

# Get the real camera locations and deployment times for use in the simulation

camera_dat <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/cameratraps/RS_2023_full_reads_all_three_cameratrap_dates.xlsx") %>% 
  mutate(time = as.numeric(difftime(time, as_datetime("1899-12-31 00:00:00"), units = "sec"))) %>% 
  filter(camera == "A", time > 60 * 10, time <= 60 * 30) # Filter out descent/retrieval frames

camera_counts <- camera_dat %>%
  group_by(Station_ID, date) %>% 
  summarize(count = sum(total), nframes = n())

camera_locations_corrected <- read_csv("intermediate/corrected_camera_stations.csv") %>% 
  select(Station_ID, Date, Longitude, Latitude)

camera_locs <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/cameratraps/RS_2023_full_reads_all_three_cameratrap_dates.xlsx", sheet = "StationData") %>% 
  filter(`Camera (A or C)` == "A") %>%
  dplyr::select(Station_ID, Date, Time = Start_Time_GMT, 
                current_dir = `Current Direction`) %>% 
  left_join(camera_counts, by = c("Station_ID", "Date" = "date")) %>% 
  left_join(camera_locations_corrected, by = c("Station_ID", "Date"))

camera_locs$current_dir <- ifelse(
  grepl("towards", tolower(camera_locs$current_dir)), 1, 
  ifelse(grepl("away", tolower(camera_locs$current_dir)), 2, 3)
)


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

# Make a matrix of observed counts
mtx_nrow <- max(table(camera_locs$Date))
ncam_vec <- numeric(3)
y_mtx <- matrix(ncol = 3, nrow = mtx_nrow)
for (i in 1:3) {
  ncam_vec[i] <- sum(camera_locs$Date == Dates[i])
  y_mtx[1:ncam_vec[i], i] <- camera_locs$count[camera_locs$Date == Dates[i]]
}

# Make a matrix of frames per camera
nframe_mtx <- y_mtx
for (i in 1:3) {
  ncam_vec[i] <- sum(camera_locs$Date == Dates[i])
  nframe_mtx[1:ncam_vec[i], i] <- camera_locs$nframes[camera_locs$Date == Dates[i]]
}

current_dir_mtx <- y_mtx 
for (i in 1:3) {
  ncam_vec[i] <- sum(camera_locs$Date == Dates[i])
  current_dir_mtx[1:ncam_vec[i], i] <- camera_locs$current_dir[camera_locs$Date == Dates[i]]
}


# Make a 3D array of camera coords  
X_array <- array(dim = c(mtx_nrow, 3, 2))
for (i in 1:3) {
  X_array[1:ncam_vec[i], i, 1] <- X$x[camera_locs$Date == Dates[i]]
  X_array[1:ncam_vec[i], i, 2] <- X$y[camera_locs$Date == Dates[i]]
}



#### Make model data lists ####

data_template <- list(
  constants = list(M = M,
                   current_dir = current_dir_mtx,
                   ncam = ncam_vec,
                   X = X_array,
                   xmin = grid_bbox[1], 
                   xmax = grid_bbox[2], 
                   ymin = grid_bbox[3], 
                   ymax = grid_bbox[4],
                   resoln = res(vps_intensity_ras)[1], 
                   hm = vps_mtx,
                   hm_nrow = nrow(vps_mtx),
                   hm_ncol = ncol(vps_mtx), 
                   nframes = nframe_mtx,
                   nROV = nrow(rov_dat),
                   ROV_hb_area = rov_dat$area * rov_dat$structured_pct
  ),
  data = list(
    y = y_mtx,
    ones = matrix(1, nrow = M, ncol = 3),
    ROV_obs = rov_dat$count
  ),
  inits = list(
    z = matrix(ncol = 3, rbinom(M*3, 1, 0.1)),
    psi = 0.1,
    lam0 = rep(0.1, 3),
    log_sigma = 3,
    s = initialize_s(n = M, t = 3,
                     xmin = grid_bbox[1],
                     xmax = grid_bbox[2],
                     ymin = grid_bbox[3],
                     ymax = grid_bbox[4],
                     resoln = 50, habitatMask = vps_mtx)
  )
)


#### nimble model code ####

my_uscr <- nimbleCode({
  
  # for (i in 1:nROV) {
    # ROV_obs[i] ~ dpois(density * ROV_hb_area[i])
  # }
  
  # Loop over all potential individuals
  for (i in 1:M) {
    # Latent state representing the inclusion prob.
    
    for (t in 1:3) {
      z[i, t] ~ dbern(psi)
      
      # Distribution of centroids
      s[i, t, 1:2] ~ dHabDistr(
        xmax = xmax,
        xmin = xmin,
        ymax = ymax,
        ymin = ymin,
        resoln = resoln,
        habitatMask = hm[1:hm_nrow, 1:hm_ncol]
      )

      # Calculate distances
      lambda[i, 1:ncam[t], t] <- GetDetectionRate(s = s[i, t, 1:2], 
                                                  X = X[1:ncam[t], t, 1:2], 
                                                  J = ncam[t], 
                                                  sigma = sigma, 
                                                  lam0 = lam0[1:3],
                                                  current_dir = current_dir[1:ncam[t], t],
                                                  z = z[i, t])
    }
  }
  
  for (t in 1:3) {
    for (j in 1:ncam[t]) {
      lambda_sum[j, t] <- sum(lambda[1:M, j, t] * nframes[j, t])
      y[j, t] ~ dnegbin(size = 1/phi_CT, prob = 1/(1 + phi_CT * lambda_sum[j, t]))
    }
    
    n[t] <- sum(z[1:M, t])
  }
  
  
  # Priors
  psi ~ dbeta(0.01, 1) # Scale prior for inclusion prob.
  for (i in 1:3) {
    lam0[i] ~ dgamma(0.1, 1)  # Detection rate per frame at centroid, depends on current dir.
  }
  phi_CT  ~ dgamma(0.05, 1)
  
  log_sigma ~ dnorm(3.435, sd = 1.138) # From telemetry
  log(sigma) <- log_sigma
})


mod <- nimbleModel(
  code = my_uscr, 
  constants = data_template$constants,
  data = data_template$data,
  inits = data_template$inits
)

cmod <- compileNimble(mod)

# We want custom samplers for everything except z, s
mcmcConf <- configureMCMC(cmod, nodes = c("z", "s", "phi_CT"))

#### Custom samplers
#use block update for lam0, psi, and var bc correlated posteriors.
mcmcConf$addSampler(target = c("lam0", "psi", "log_sigma"),
                    type = 'AF_slice', control = list(adaptive=TRUE), silent = TRUE)

mcmcConf$setMonitors(c("lam0", "psi", "n", "log_sigma", "sigma", "phi_CT"))


mcmc <- buildMCMC(mcmcConf)
cmcmc <- compileNimble(mcmc)


# samples <- runMCMC(cmcmc, niter = 1000, nburnin = 0, nchains = 2,
samples <- runMCMC(cmcmc, niter = 50000, nburnin = 10000, nchains = 2,
                   thin = 2, samplesAsCodaMCMC = TRUE, 
                   inits = data_template$inits)
# samples <- runMCMC(cmcmc, niter = 2000, nburnin = 0, nchains = 1,
#                    thin = 1, samplesAsCodaMCMC = TRUE)

summary <- MCMCvis::MCMCsummary(samples)
summary$param <- rownames(summary)

rownames(summary) <- NULL

saveRDS(list(summary = summary,
             samples = samples),
        paste0("uSCR_real/joint_masked_VPSsurface_NB.RDS"))




