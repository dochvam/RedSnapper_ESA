library(tidyverse)
library(readxl)
library(terra)
library(parallel)
library(nimble)
library(coda)
library(MCMCvis)

set.seed(4077)

# Source file: Ben Augustine's helper fns
source("uSCR_binom_Augustine/augustine_helper.R")

# Source file: Ben Goldstein's helper fns
source("uSCR_binom_Augustine/other_helper.R")

M <- 3000
integration_type <- "Full"
spatial_type <- "HB_map"
source('pipeline_FL/prep_data_FL.R', local = TRUE)


model_code <- nimbleCode({
  lambda.N ~ dunif(0, 50000) #expected abundance
  N ~ dpois(lambda.N) #realized abundance
  spatial_beta ~ dnorm(1, sd = 1)
  theta ~ dunif(0, 10)
  
  phi[1:hm_nrow, 1:hm_ncol] <- exp(spatial_beta * log(hm[1:hm_nrow, 1:hm_ncol])) # hm is log-scale covariate
  
  
  if (integration_type == "Camera_ROV" || integration_type == "Full") {
    for (i in 1:nROV) {
      # rbs and rbe are ROV buffer start/end indexes for ROV i
      pctFishInROVbuffer[i] <- calcPctFishInROVbuffer(phi = phi[1:hm_nrow, 1:hm_ncol], 
                                                      weights = rb_weights[rbs[i]:rbe[i]], 
                                                      rov_cell_xvec = rov_cell_xvec[rbs[i]:rbe[i]],
                                                      rov_cell_yvec = rov_cell_yvec[rbs[i]:rbe[i]],
                                                      n = rbe[i] - rbs[i] + 1)
      
      if (distr == "Pois") {
        ROV_obs[i] ~ dpois(pctFishInROVbuffer[i] * lambda.N)
      } else if (distr == "NB") {
        ROV_obs[i] ~ dnegbin(size = 1/theta, prob = 1/(1 + theta * pctFishInROVbuffer[i] * lambda.N))
      }
    }
  }
}) #model


mod <- nimbleModel(model_code, 
                   constants = list(
                     distr = "Pois",
                     hm_nrow = constants$hm_nrow,
                     hm_ncol = constants$hm_ncol,
                     hm = constants$hm,
                     nROV = constants$nROV,
                     rov_cell_xvec = constants$rov_cell_xvec,
                     rov_cell_yvec = constants$rov_cell_yvec,
                     rb_weights = constants$rb_weights,
                     rbe = constants$rbe,
                     rbs = constants$rbs
                   ),
                   inits = list(
                     spatial_beta = 1, lambda.N = 1000, N = 1000, theta = 0.1
                   ),
                   data = list(ROV_obs = Nimdata$ROV_obs)
)
conf <- configureMCMC(mod)

mcmc <- buildMCMC(conf)

clist <- compileNimble(mod, mcmc)

mcmc_samples <- runMCMC(clist$mcmc, niter = 40000, nburnin = 10000, nchains = 2,
                        thin = 1, samplesAsCodaMCMC = TRUE)

summary <- MCMCvis::MCMCsummary(mcmc_samples)

saveRDS(list(
  samples = mcmc_samples,
  summary = summary
), paste0("pipeline_FL/FL_results/ROV_only_", spatial_type, ".RDS"))
