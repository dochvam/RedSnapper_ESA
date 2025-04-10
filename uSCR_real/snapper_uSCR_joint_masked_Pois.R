### In this script I modify the "simple" model to include:
#' - Auxilliary information from ROVs, an independent source of data on density
#' - A habitat mask, only allowing fish centroids to occur on hardbottom

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


dMyHabMask <- nimbleFunction(run = function(x = double(0),
                                            s = double(1),
                                            xmax = double(0),
                                            xmin = double(0),
                                            ymax = double(0),
                                            ymin = double(0),
                                            resoln = double(0),
                                            habitatMask = double(2),
                                            log = logical(0)
){
  if (s[1] < xmin) return(-Inf)
  if (s[1] > xmax) return(-Inf)
  if (s[2] < ymin) return(-Inf)
  if (s[2] > ymax) return(-Inf)
  
  # Convert s to row/col
  row <- trunc((s[1] - xmin) / resoln) + 1
  col <- trunc((s[2] - ymin) / resoln) + 1
  
  test <- 1 - (habitatMask[row, col] == 0)
  
  if (log) return(log(test))
  else return(test)
  returnType(double(0))
})

rMyHabMask <- nimbleFunction(run = function(n = double(0),
                                            s = double(1),
                                            xmax = double(0),
                                            xmin = double(0),
                                            ymax = double(0),
                                            ymin = double(0),
                                            resoln = double(0),
                                            habitatMask = double(2)
){
  returnType(double(0))
  
  if (s[1] < xmin) return(0)
  if (s[1] > xmax) return(0)
  if (s[2] < ymin) return(0)
  if (s[2] > ymax) return(0)
  
  # Convert s to row/col
  row <- trunc((s[1] - xmin) / resoln) + 1
  col <- trunc((s[2] - ymin) / resoln) + 1
  
  if (habitatMask[row, col] == 0)  {
    return(0)
  } else {
    return(1)
  }
})

initialize_s <- function(M, habitatMask, resoln, xmin, ymin) {
  valid_indexes <- which(habitatMask == 1, arr.ind = TRUE)
  
  s_init <- array(NA, dim = c(M, 3, 2))
  for (i in 1:M) {
    for (t in 1:3) {
      ind <- valid_indexes[sample(1:nrow(valid_indexes), size = 1), ]
      s_init[i, t, ] <- c(
        xmin + resoln * ind[1] - resoln/2,
        ymin + resoln * ind[2] - resoln/2
      )
    }
  }
  s_init
}

#### Construct the habitat mask ####

hb_raw <- rast("Data/Chicken_Rock_Map/ChickenRock_Classification.tif")

terra::values(hb_raw) <- ifelse(terra::values(hb_raw) == 4, 1, 0)
terra::values(hb_raw)[is.na(terra::values(hb_raw))] <- 0

hb_mask <- hb_raw %>% 
  aggregate(30, fun = "max")
terra::values(hb_mask)[is.nan(terra::values(hb_mask))] <- 0


hb_mtx <- t(as.matrix(hb_mask, wide = T))
hb_mtx <- hb_mtx[, ncol(hb_mtx):1]
hb_bbox <- ext(hb_mask)
hb_bbox <- as.numeric(c(hb_bbox[1], hb_bbox[2], hb_bbox[3], hb_bbox[4]))

x_offset <- mean(hb_bbox[1:2])
y_offset <- mean(hb_bbox[3:4])

hb_bbox[1:2] <- hb_bbox[1:2] - x_offset
hb_bbox[3:4] <- hb_bbox[3:4] - y_offset


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
M <- 1000

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
  project(crs(hb_mask)) %>% 
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



#### Make model data lists

data_template <- list(
  constants = list(M = M,
                   current_dir = current_dir_mtx,
                   ncam = ncam_vec,
                   X = X_array,
                   xmin = hb_bbox[1], 
                   xmax = hb_bbox[2], 
                   ymin = hb_bbox[3], 
                   ymax = hb_bbox[4],
                   resoln = res(hb_mask)[1],
                   hb = hb_mtx,
                   hb_area = sum(hb_mtx) * res(hb_mask)[1]^2,
                   nframes = nframe_mtx,
                   nxcell = nrow(hb_mtx), 
                   nycell = ncol(hb_mtx),
                   nROV = nrow(rov_dat),
                   ROV_hb_area = rov_dat$area * rov_dat$structured_pct),
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
    s = initialize_s(M = M, habitatMask = hb_mtx, 
                     xmin = hb_bbox[1], ymin = hb_bbox[3],
                     resoln = res(hb_mask)[1])
  )
)


#### nimble model code ####

my_uscr <- nimbleCode({
  
  for (i in 1:nROV) {
    ROV_obs[i] ~ dpois(density * ROV_hb_area[i])
  }
  
  # Loop over all potential individuals
  for (i in 1:M) {
    # Latent state representing the inclusion prob.
    
    for (t in 1:3) {
      z[i, t] ~ dbern(psi)
      
      # Distribution of centroids
      s[i, t, 1] ~ dunif(xmin, xmax)
      s[i, t, 2] ~ dunif(ymin, ymax)
      # Ones trick to reject s outside the habitat mask
      ones[i, t] ~ dMyHabMask(
        s = s[i, t, 1:2], xmax = xmax, xmin = xmin, ymax = ymax, ymin = ymin,
        resoln = resoln, habitatMask = hb[1:nxcell, 1:nycell]
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
      y[j, t] ~ dpois(lambda_sum[j, t])
    }
    
    n[t] <- sum(z[1:M, t])
  }
  
  
  # Priors
  psi ~ dbeta(0.01, 1) # Scale prior for inclusion prob.
  for (i in 1:3) {
    lam0[i] ~ dgamma(0.1, 1)  # Detection rate per frame at centroid, depends on current dir.
  }
  density <- psi * M / hb_area
  
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
mcmcConf <- configureMCMC(cmod, nodes = c("z", "s"))

#### Custom samplers
#use block update for lam0, psi, and var bc correlated posteriors.
mcmcConf$addSampler(target = c("lam0", "psi", "log_sigma"),
                    type = 'AF_slice', control = list(adaptive=TRUE), silent = TRUE)

mcmcConf$setMonitors(c("lam0", "psi", "n", "log_sigma", "sigma", "density"))


mcmc <- buildMCMC(mcmcConf)
cmcmc <- compileNimble(mcmc)


# samples <- runMCMC(cmcmc, niter = 1000, nburnin = 0, nchains = 2,
samples <- runMCMC(cmcmc, niter = 20000, nburnin = 5000, nchains = 2,
                   thin = 5, samplesAsCodaMCMC = TRUE, 
                   inits = data_template$inits)
# samples <- runMCMC(cmcmc, niter = 2000, nburnin = 0, nchains = 1,
#                    thin = 1, samplesAsCodaMCMC = TRUE)

summary <- MCMCvis::MCMCsummary(samples)
summary$param <- rownames(summary)

rownames(summary) <- NULL

saveRDS(list(summary = summary,
             samples = samples),
        paste0("uSCR_real/joint_masked_posterior_Pois.RDS"))




