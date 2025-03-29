library(nimble)
library(tidyverse)
library(readxl)
library(terra)

dZIP <- nimbleFunction(
  run = function(x = double(), lambda = double(),
                 zeroProb = double(), log = logical(0, default = 0)) {
    returnType(double())
    ## For use with AD, we cannot use an `if` statement to handle the mixture.
    prob <- zeroProb * dbinom(x, size = 1, prob = 0) + (1 - zeroProb) * dpois(x, lambda)
    if (log) return(log(prob))
    return(prob)
  },
  buildDerivs = 'run'   # Needed when used with AD-based algorithms.
)
rZIP <- nimbleFunction(
  run = function(n = integer(), lambda = double(), zeroProb = double()) {
    returnType(double())
    isStructuralZero <- rbinom(1, prob = zeroProb, size = 1)
    if (isStructuralZero) return(0)
    return(rpois(1, lambda))
  })

#### Load the camera data ####

camera_dat <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/cameratraps/RS_2023_full_reads_all_three_cameratrap_dates.xlsx") %>% 
  mutate(time = as.numeric(difftime(time, as_datetime("1899-12-31 00:00:00"), units = "sec"))) %>% 
  filter(camera == "A", time > 60 * 10 & time <= 60 * 30) # Filter out descent/retrieval frames


camera_locations_corrected <- read_csv("intermediate/corrected_camera_stations.csv") %>% 
  select(Station_ID, Date, Longitude, Latitude)

stations <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/cameratraps/RS_2023_full_reads_all_three_cameratrap_dates.xlsx", sheet = "StationData") %>% 
  mutate(Start_Time_GMT = as.numeric(difftime(Start_Time_GMT, as_datetime("1899-12-31 00:00:00"), units = "hours"))) %>% 
  filter(`Camera (A or C)` == "A") %>% 
  select(-Start_Longitude, -Start_Latitude) %>% 
  left_join(camera_locations_corrected)

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
         hb_area = area * structured_pct)


#### Rearrange data for model ####
dethists <- camera_dat %>% 
  dplyr::select(date, Station_ID, total) %>% 
  group_by(date, Station_ID) %>% 
  mutate(frame = paste0("F", row_number())) %>% 
  ungroup() %>% 
  pivot_wider(values_from = total, names_from = frame) %>% 
  left_join(distinct(stations, date = Date, Station_ID, 
                     Latitude, Longitude,
                     current_dir = `Current Direction`), by = c("date", "Station_ID"))

dethist_mtx <- as.matrix(dethists[, grep("^F", colnames(dethists))])
y_cam <- unlist(apply(dethist_mtx, 1, sum, na.rm = T))
nframes <- rowSums(!is.na(dethist_mtx))

current_dir <- ifelse(
  grepl("towards", tolower(dethists$current_dir)), 1, 
  ifelse(grepl("away", tolower(dethists$current_dir)), 2, 3)
)

#### Get distance to HB ####

hardbottom <- rast("Data/Chicken_Rock_Map/ChickenRock_Classification.tif")
terra::values(hardbottom) <- as.numeric(terra::values(hardbottom) == 4)


rov_pts <- vect(rov_dat, geom = c("Longitude", "Latitude"),
                crs = "+proj=longlat") %>% 
  project(crs(hardbottom))
rov_buffers <- buffer(rov_pts, 100)
pct_hb_ROV <- numeric(nrow(rov_buffers))
for (i in 1:nrow(rov_buffers)) {
  pct_hb_ROV[i] <- mean(terra::extract(hardbottom, rov_buffers[i, ])[, 2], na.rm = T)
}

cam_pts <- vect(dethists, geom = c("Longitude", "Latitude"),
                           crs = "+proj=longlat") %>% 
  project(crs(hardbottom))
cam_buffers <- buffer(cam_pts, 100)
pct_hb_cam <- numeric(nrow(cam_buffers))
for (i in 1:nrow(cam_buffers)) {
  pct_hb_cam[i] <- mean(terra::extract(hardbottom, cam_buffers[i, ])[, 2], na.rm = T)
}



#### Define a simple model ####
rov_cam_code <- nimbleCode({
  for (i in 1:ncamera) {
    log(lambda_cam[i]) <- log(lambda0) + beta_hb * pct_hb_cam[i]
    N[i] ~ dpois(lambda_cam[i] * theta[current_dir[i]])
    for (j in 1:nframes[i]) {
      y_cam[i, j] ~ dbinom(size = N[i], prob = p_nmix)
    }
  }
  
  for (i in 1:nROV) {
    log(lambda_ROV[i]) <- log(lambda0) + beta_hb * pct_hb_ROV[i]
    
    y_ROV[i] ~ dpois(rov_area[i] * lambda_ROV[i])
  }
  
  for (i in 1:3) {
    theta[i] ~ dunif(0, 1000000)
  }
  
  lambda0 ~ dgamma(1, 0.5)
  p_nmix ~ dunif(0, 1)
  beta_hb ~ dnorm(0, sd = 10)
})

rov_cam_mod_nmix <- nimbleModel(
  code = rov_cam_code,
  data = list(
    y_cam = dethist_mtx,
    y_ROV = rov_dat$count
  ),
  constants = list(
    current_dir = current_dir,
    rov_area = rov_dat$area,
    ncamera = nrow(dethist_mtx),
    nROV = nrow(rov_dat),
    pct_hb_ROV = pct_hb_ROV,
    pct_hb_cam = pct_hb_cam,
    nframes = nframes
  ),
  inits = list(lambda0 = 0.1, 
               beta_hb = 0, p_nmix = 0.1,
               N = rep(max(dethist_mtx) + 1, nrow(dethist_mtx)),
               theta = rep(100, 3))
)

conf <- configureMCMC(rov_cam_mod_nmix)

mcmc <- buildMCMC(conf)

complist <- compileNimble(rov_cam_mod_nmix, mcmc)

samples <- runMCMC(complist$mcmc, niter = 10000, nburnin = 2000,
                   nchains = 3, thin = 5, samplesAsCodaMCMC = TRUE)

summary <- MCMCvis::MCMCsummary(samples)
summary$param <- rownames(summary)
summary$distr <- "Nmix_Poi"
rownames(summary) <- NULL

write_csv(summary, "offset_counts/Pois_Nmix_offset_results.csv")

summary %>% 
  left_join(data.frame(
    param = paste0("theta[", 1:3, "]"),
    dir = c("Towards", "Away", "Perpendicular")
  )) %>% 
  filter(!is.na(dir)) %>% 
  ggplot() +
  geom_pointrange(aes(dir, mean, ymin = `2.5%`, ymax = `97.5%`)) +
  ylab("Effective sampling area") + xlab("Current dir.") +
  coord_flip() +
  theme_minimal()


#### Repeat with a BBP N-mixture ####

rov_cam_code <- nimbleCode({
  for (i in 1:ncamera) {
    log(lambda_cam[i]) <- log(lambda0) + beta_hb * pct_hb_cam[i]
    y_cam[i, 1:nframes[i]] ~ dNmixture_BBP_s(
      lambda = lambda_cam[i] * theta[current_dir[i]], 
      # theta = phi, 
      s = s, prob = p_nmix, Nmin = 0, Nmax = 1000, 
      len = nframes[i]
    )
  }
  
  for (i in 1:nROV) {
    log(lambda_ROV[i]) <- log(lambda0) + beta_hb * pct_hb_ROV[i]
    
    y_ROV[i] ~ dpois(rov_area[i] * lambda_ROV[i])
  }
  
  for (i in 1:3) {
    theta[i] ~ dunif(0, 1000000)
  }
  
  # phi ~ dgamma(1, 0.5)
  s ~ dgamma(1, 0.5)
  lambda0 ~ dgamma(1, 0.5)
  p_nmix ~ dunif(0, 1)
  beta_hb ~ dnorm(0, sd = 10)
})

rov_cam_mod_nmix <- nimbleModel(
  code = rov_cam_code,
  data = list(
    y_cam = dethist_mtx,
    y_ROV = rov_dat$count
  ),
  constants = list(
    current_dir = current_dir,
    rov_area = rov_dat$area,
    ncamera = nrow(dethist_mtx),
    nROV = nrow(rov_dat),
    pct_hb_ROV = pct_hb_ROV,
    pct_hb_cam = pct_hb_cam,
    nframes = nframes
  ),
  inits = list(lambda0 = 0.1, 
               s = 0.1,
               beta_hb = 0, p_nmix = 0.1,
               theta = rep(100, 3))
)

conf <- configureMCMC(rov_cam_mod_nmix)

mcmc <- buildMCMC(conf)

complist <- compileNimble(rov_cam_mod_nmix, mcmc)

samples <- runMCMC(complist$mcmc, niter = 4000, nburnin = 2000,
                   nchains = 2, thin = 1, samplesAsCodaMCMC = TRUE)

summary <- MCMCvis::MCMCsummary(samples)
summary$param <- rownames(summary)
summary$distr <- "Nmix_BBP"
rownames(summary) <- NULL

write_csv(summary, "offset_counts/BBP_Nmix_offset_results.csv")

bind_rows(
    read_csv("offset_counts/BBP_Nmix_offset_results.csv"),
    read_csv("offset_counts/Pois_model_results.csv"),
    read_csv("offset_counts/ZIP_model_results.csv"),
    read_csv("offset_counts/NB_model_results.csv"),
    read_csv("offset_counts/Pois_Nmix_offset_results.csv")
  ) %>% 
  left_join(data.frame(
    param = paste0("theta[", 1:3, "]"),
    dir = c("Towards", "Away", "Perpendicular")
  )) %>% 
  filter(!is.na(dir)) %>% 
  ggplot() +
  geom_pointrange(aes(dir, mean, col = distr, ymin = `2.5%`, ymax = `97.5%`),
                  position = position_dodge(width = 0.2)) +
  ylab("Effective sampling area") + xlab("Current dir.") +
  coord_flip() +
  theme_minimal()



