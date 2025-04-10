library(nimbleEcology)
library(tidyverse)
library(readxl)
library(terra)


#### Construct the habitat mask ####

hb_raw <- rast("Data/Chicken_Rock_Map/ChickenRock_Classification.tif")

terra::values(hb_raw) <- ifelse(terra::values(hb_raw) == 4, 1, 0)
terra::values(hb_raw)[is.na(terra::values(hb_raw))] <- 0

hb_mask <- hb_raw %>% 
  aggregate(30, fun = "max")
terra::values(hb_mask)[is.nan(terra::values(hb_mask))] <- 0


hb_mtx <- t(as.matrix(hb_mask, wide = T))
hb_mtx <- hb_mtx[, ncol(hb_mtx):1]



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

#### What does a regression of ROV look like? ####

rov_only_code <- nimbleCode({
  for (i in 1:nROV) {
    ROV_obs[i] ~ dnegbin(size = 1/phi, prob = 1 / (1 + phi * density * ROV_hb_area[i]))
  }
  
  psi ~ dbeta(0.01, 1) # Use the same prior
  phi ~ dgamma(0.5, 1)
  density <- psi * M / hb_area
})

mod <- nimbleModel(
  code = rov_only_code, 
  constants = list(
    M = 5000,
    nROV = nrow(rov_dat),
    ROV_hb_area = rov_dat$area * rov_dat$structured_pct,
    hb_area = sum(hb_mtx) * res(hb_mask)[1]^2
  ),
  data = list(
    ROV_obs = rov_dat$count
  ),
  inits = list(
    psi = 0.1,
    phi = 0.1
  )
)

samples_ROVonly <- nimbleMCMC(mod, niter = 10000, nburnin = 5000, nchains = 2,
                              samplesAsCodaMCMC = TRUE,
                              monitors = c("psi", "density", "phi"))
summary_ROVonly <- MCMCvis::MCMCsummary(samples_ROVonly)

saveRDS(list(summary = summary_ROVonly,
             samples = samples_ROVonly),
        paste0("uSCR_real/ROV_only_posterior.RDS"))
