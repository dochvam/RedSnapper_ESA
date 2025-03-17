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
  filter(camera == "A", time >= 60 * 10 & time < 60 * 30) # Filter out descent/retrieval frames

stations <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/cameratraps/RS_2023_full_reads_all_three_cameratrap_dates.xlsx", sheet = "StationData") %>% 
  mutate(Start_Time_GMT = as.numeric(difftime(Start_Time_GMT, as_datetime("1899-12-31 00:00:00"), units = "hours"))) %>% 
  filter(`Camera (A or C)` == "A") %>% 
  mutate(Start_Longitude = as.numeric(Start_Longitude))

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
  mutate(density = count/area)


#### Rearrange data for model ####
dethists <- camera_dat %>% 
  dplyr::select(date, Station_ID, total) %>% 
  group_by(date, Station_ID) %>% 
  mutate(frame = paste0("F", row_number())) %>% 
  ungroup() %>% 
  pivot_wider(values_from = total, names_from = frame) %>% 
  left_join(distinct(stations, date = Date, Station_ID, 
                     Latitude = Start_Latitude, Longitude = Start_Longitude,
                     current_dir = `Current Direction`), by = c("date", "Station_ID"))

dethist_mtx <- as.matrix(dethists[, grep("^F", colnames(dethists))])
y_cam <- unlist(apply(dethist_mtx, 1, mean, na.rm = T))

current_dir <- ifelse(
  grepl("towards", tolower(dethists$current_dir)), 1, 
  ifelse(grepl("away", tolower(dethists$current_dir)), 2, 3)
)

#### Get distance to HB ####

hardbottom <- rast("Data/Chicken_Rock_Map/ChickenRock_Classification.tif")
values(hardbottom) <- as.numeric(values(hardbottom) == 4)


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
    if (distr == "NB") {
      y_cam[i] ~ dnegbin(size = 1/phi1, prob = 1/(1 + phi1 * theta[current_dir[i]] * lambda_cam[i]))
      y_cam_ppc[i] ~ dnegbin(size = 1/phi1, prob = 1/(1 + phi1 * theta[current_dir[i]] * lambda_cam[i]))
    } else if (distr == "ZIP") {
      y_cam[i] ~ dZIP(lambda = theta[current_dir[i]] * lambda_cam[i], zeroProb = zipCam)
      y_cam_ppc[i] ~ dZIP(lambda = theta[current_dir[i]] * lambda_cam[i], zeroProb = zipCam)
    } else if (distr == "Pois") {
      y_cam[i] ~ dpois(theta[current_dir[i]] * lambda_cam[i])
      y_cam_ppc[i] ~ dpois(theta[current_dir[i]] * lambda_cam[i])
    }
  }
  cam_var <- var(y_cam_ppc[1:ncamera])
  
  for (i in 1:nROV) {
    log(lambda_ROV[i]) <- log(lambda0) + beta_hb * pct_hb_ROV[i]
    if (distr == "NB") {
      y_ROV[i] ~ dnegbin(size = 1/phi2, prob = 1/(1 + phi2 * rov_area[i] * lambda_ROV[i]))
    } else if (distr == "ZIP") {
      y_ROV[i] ~ dZIP(lambda = rov_area[i] * lambda_ROV[i], zeroProb = zipROV)
    } else if (distr == "Pois") {
      y_ROV[i] ~ dpois(rov_area[i] * lambda_ROV[i])
    }
  }
  
  for (i in 1:3) {
    theta[i] ~ dunif(0, 10000)
  }
  
  lambda0 ~ dgamma(1, 0.5)
  beta_hb ~ dnorm(0, sd = 10)
  zipROV ~ dunif(0, 1)
  zipCam ~ dunif(0, 1)
  phi1 ~ dgamma(1, 0.5)
  phi2 ~ dgamma(1, 0.5)
})

rov_cam_mod_simple <- nimbleModel(
  code = rov_cam_code,
  data = list(
    y_cam = y_cam,
    y_ROV = rov_dat$count
  ),
  constants = list(
    current_dir = current_dir,
    rov_area = rov_dat$area,
    ncamera = length(y_cam),
    nROV = nrow(rov_dat),
    pct_hb_ROV = pct_hb_ROV,
    pct_hb_cam = pct_hb_cam,
    distr = "Pois"
  ),
  inits = list(lambda0 = 0.1, 
               beta_hb = 0,
               theta = rep(100, 3), 
               phi1 = 10, phi2 = 10)
)

conf <- configureMCMC(rov_cam_mod_simple)
conf$addMonitors("cam_var")

mcmc <- buildMCMC(conf)

complist <- compileNimble(rov_cam_mod_simple, mcmc)

samples <- runMCMC(complist$mcmc, niter = 10000, nburnin = 2000,
                   nchains = 3, thin = 5, samplesAsCodaMCMC = TRUE)

summary <- MCMCvis::MCMCsummary(samples)
summary$param <- rownames(summary)
summary$distr <- "Pois"
rownames(summary) <- NULL

write_csv(summary, "intermediate/simple_model_results.csv")

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

#### PPC: are the data overdispersed for a Poisson? #### 
cam_var_samples <- as.numeric(unlist(samples[, "cam_var"]))
mean(var(y_cam) > cam_var_samples)

{
  hist(cam_var_samples)
  abline(v = var(y_cam))
}

#### Run it again, this time using a NB ####

rov_cam_mod_simple <- nimbleModel(
  code = rov_cam_code,
  data = list(
    y_cam = y_cam,
    y_ROV = rov_dat$count
  ),
  constants = list(
    current_dir = current_dir,
    rov_area = rov_dat$area,
    ncamera = length(y_cam),
    nROV = nrow(rov_dat),
    pct_hb_ROV = pct_hb_ROV,
    pct_hb_cam = pct_hb_cam,
    distr = "NB"
  ),
  inits = list(lambda0 = 0.1, 
               beta_hb = 0,
               theta = rep(100, 3), 
               phi1 = 10, phi2 = 10)
)

conf <- configureMCMC(rov_cam_mod_simple)
conf$addMonitors("cam_var")
mcmc <- buildMCMC(conf)

complist <- compileNimble(rov_cam_mod_simple, mcmc)

samplesNB <- runMCMC(complist$mcmc, niter = 10000, nburnin = 2000,
                   nchains = 3, thin = 5, samplesAsCodaMCMC = TRUE)

summaryNB <- MCMCvis::MCMCsummary(samplesNB)
summaryNB$param <- rownames(summaryNB)
summaryNB$distr <- "NB"
rownames(summaryNB) <- NULL

write_csv(bind_rows(summary, summaryNB), "intermediate/simple_model_results.csv")

bind_rows(summary, summaryNB) %>% 
  left_join(data.frame(
    param = paste0("theta[", 1:3, "]"),
    dir = c("Towards", "Away", "Perpendicular")
  )) %>% 
  filter(!is.na(dir)) %>% 
  ggplot() +
  geom_pointrange(aes(dir, mean, ymin = `2.5%`, ymax = `97.5%`,
                      col = distr), position = position_dodge(width = 0.1)) +
  ylab("Effective sampling area") + xlab("Current dir.") +
  coord_flip() +
  theme_minimal()

#### PPC: are the data overdispersed for a NB? #### 
cam_var_samples <- as.numeric(unlist(samplesNB[, "cam_var"]))
mean(var(y_cam) > cam_var_samples)

{
  hist(cam_var_samples)
  abline(v = var(y_cam))
}

#### Run it again, this time using a ZIP ####

rov_cam_mod_simple <- nimbleModel(
  code = rov_cam_code,
  data = list(
    y_cam = y_cam,
    y_ROV = rov_dat$count
  ),
  constants = list(
    current_dir = current_dir,
    rov_area = rov_dat$area,
    ncamera = length(y_cam),
    nROV = nrow(rov_dat),
    pct_hb_ROV = pct_hb_ROV,
    pct_hb_cam = pct_hb_cam,
    distr = "ZIP"
  ),
  inits = list(lambda0 = 0.1, 
               beta_hb = 0,
               theta = rep(100, 3), 
               phi1 = 10, phi2 = 10)
)

conf <- configureMCMC(rov_cam_mod_simple)
conf$addMonitors("cam_var")
mcmc <- buildMCMC(conf)

complist <- compileNimble(rov_cam_mod_simple, mcmc)

samplesZIP <- runMCMC(complist$mcmc, niter = 10000, nburnin = 2000,
                   nchains = 3, thin = 5, samplesAsCodaMCMC = TRUE)

summaryZIP <- MCMCvis::MCMCsummary(samplesZIP)
summaryZIP$param <- rownames(summaryZIP)
summaryZIP$distr <- "ZIP"
rownames(summaryZIP) <- NULL

write_csv(bind_rows(summary, summaryNB, summaryZIP), "intermediate/simple_model_results.csv")

bind_rows(summary, summaryNB, summaryZIP) %>% 
  left_join(data.frame(
    param = paste0("theta[", 1:3, "]"),
    dir = c("Towards", "Away", "Perpendicular")
  )) %>% 
  filter(!is.na(dir)) %>% 
  ggplot() +
  geom_pointrange(aes(dir, mean, ymin = `2.5%`, ymax = `97.5%`,
                      col = distr), position = position_dodge(width = 0.1)) +
  ylab("Effective sampling area") + xlab("Current dir.") +
  coord_flip() +
  theme_minimal()

#### PPC: are the data overdispersed for a ZIP? #### 
cam_var_samples <- as.numeric(unlist(samplesZIP[, "cam_var"]))
mean(var(y_cam) > cam_var_samples)

{
  hist(cam_var_samples)
  abline(v = var(y_cam))
}

