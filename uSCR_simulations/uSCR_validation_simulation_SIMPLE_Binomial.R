library(tidyverse)
library(readxl)
library(nimble)

##########################################################
#' Do we get anything useful out of our SCR model at the
#' scale we're working at?
#' 
#' Update "SIMPLE" to use binomial formulation for fish
##########################################################
source("uSCR_simulations/SCR_sim_fn.R")

#### Parallel execution ####
library(parallel)
cl <- makeCluster(4)
captured <- clusterEvalQ(cl, source("uSCR_simulations/SCR_sim_fn.R"))

result <- parLapply(cl = cl, X = 600 + 1:4, fun = run_one_uSCR_simulation_binomial,
                    niter = 5000, nburnin = 0, nchains = 2, thin = 1,
                    M = 50, sampler_spec = "RJMCMC", prefix = "sim_validation_RJMCMC_")

summary <- list.files("intermediate/sim/", pattern = "sim_validation_RJMCMC_simplebinomial_", full.names = T) %>% 
  lapply(function(x) readRDS(x)$summary) %>% 
  bind_rows()

summary %>% 
  filter(param == "n") %>% 
  ggplot() + 
  geom_point(aes(iter, true_N/3), col = "red", size = 2) +
  geom_pointrange(aes(iter, mean, ymin = `2.5%`, ymax = `97.5%`)) +
  ylab("Estimated population") + xlab("Simulation iter.") +
  theme_minimal() +
  coord_flip()
summary %>% 
  filter(param == "log_sigma") %>% 
  ggplot() + 
  # geom_point(aes(iter, true_N), col = "red", size = 2) +
  geom_pointrange(aes(iter, mean, ymin = `2.5%`, ymax = `97.5%`)) +
  ylab("Estimate of log_sigma") + xlab("Simulation iter.") +
  theme_minimal() +
  coord_flip()
summary %>% 
  filter(param == "p0") %>% 
  ggplot() + 
  # geom_point(aes(iter, true_N), col = "red", size = 2) +
  geom_pointrange(aes(iter, mean, ymin = `2.5%`, ymax = `97.5%`)) +
  ylab("Estimate of p0") + xlab("Simulation iter.") +
  theme_minimal() +
  coord_flip()


# Check mixing of s, z


temp <- readRDS("intermediate/sim/sim_validation_RJMCMC_simplebinomial_uSCR_301.RDS")

data.frame(x = as.numeric(unlist(temp$samples[[1]][,"s[109, 1, 1]"])),
     y = as.numeric(unlist(temp$samples[[1]][,"s[109, 1, 2]"])),
     z = as.numeric(unlist(temp$samples[[1]][,"z[109, 1]"]))) %>% 
  ggplot() + 
  geom_point(aes(x, y, col = z))

