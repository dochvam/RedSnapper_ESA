library(tidyverse)
library(nimble)
library(terra)
library(readxl)
library(measurements)

iter <- 10
M <- 900
integration_type <- "Full"
source("uSCR_binom_Augustine/other_helper.R")
source("pipeline_NC/prep_data_NC.R")


# Load results from camera-only model
# cam_summary <- readRDS("pipeline_NC/NC_results/uSCR_real_Augustine_Binom_LargerBuffer_100_Camera_Telemetry.RDS")$summary
# full_summary <- readRDS("pipeline_NC/NC_results/uSCR_real_Augustine_Binom_LargerBuffer_100_Camera_ROV.RDS")$summary

results_files_full <- list.files("pipeline_NC/NC_results/", 
                                 pattern = "uSCR_real_Augustine_Binom_20minSigma_*.*Full.RDS",
                                 full.names = TRUE)

samples_full <- lapply(results_files_full, function(x) readRDS(x)$samples) %>% 
  as.mcmc.list()

full_summary <- MCMCsummary(samples_full) %>% 
  mutate(integration_type = "Full")
full_summary$param <- rownames(full_summary)
rownames(full_summary) <- NULL


results_files_cam <- list.files("pipeline_NC/NC_results/", 
                                pattern = "uSCR_real_Augustine_Binom_20minSigma_*.*Camera_only.RDS",
                                full.names = TRUE)

samples_cam <- lapply(results_files_cam, function(x) readRDS(x)$samples) %>% 
  as.mcmc.list()

cam_summary <- MCMCsummary(samples_cam) %>% 
  mutate(integration_type = "Camera_only")
cam_summary$param <- rownames(cam_summary)
rownames(cam_summary) <- NULL


# Plot the point estimate of fish per phi
spatial_beta <- 0.5
N <- 3703.2239918 / 3

phi_unscaled <- exp(spatial_beta * log(vps_mtx))
phi <- phi_unscaled / sum(phi_unscaled)

# Expected num. fish per cell
fish_percell <- as.numeric(phi * N)
density_percell <- fish_percell / (resoln^2)
hist(density_percell)

# Get estimates of density from ROV
rov_dens <- rov_dat$count / rov_dat$area

data.frame(
  density = c(density_percell, rov_dens),
  type = c(rep("USCR estimate", length(density_percell)), rep("ROV naive", length(rov_dens)))
) %>% 
  ggplot() +
  geom_violin(aes(type, log(density + 0.0000001), group = type, col = type)) +
  geom_point(aes(type, log(density + 0.0000001), col = type),
             position = position_jitter(width = 0.2)) +
  coord_flip()


#### What counts do we *expect* at the ROVs, given camera model? ####

expectedCount <- numeric(nrow(rov_dat))

for (i in 1:nrow(rov_dat)) {
  expectedCount[i] <- N * calcPctFishInROVbuffer(phi = phi, 
                                             weights = constants$rb_weights[constants$rbs[i]:constants$rbe[i]], 
                                             rov_cell_xvec = constants$rov_cell_xvec[constants$rbs[i]:constants$rbe[i]],
                                             rov_cell_yvec = constants$rov_cell_yvec[constants$rbs[i]:constants$rbe[i]],
                                             n = constants$rbe[i] - constants$rbs[i] + 1)
}

plot(expectedCount, rov_dat$count)
