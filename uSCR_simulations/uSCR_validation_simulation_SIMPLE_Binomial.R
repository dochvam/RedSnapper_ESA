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
cl <- makeCluster(8)
captured <- clusterEvalQ(cl, source("uSCR_simulations/SCR_sim_fn.R"))

result <- parLapply(cl = cl, X = 200 + 1:8, fun = run_one_uSCR_simulation_binomial,
                    niter = 500, nburnin = 0, nchains = 2, thin = 1,
                    M = 600, sampler_spec = "RW_block_1", prefix = "sim_validation_")

summary <- list.files("intermediate/sim/", pattern = "sim_validation_simplebinomial_", full.names = T) %>% 
  lapply(function(x) readRDS(x)$summary) %>% 
  bind_rows()

summary %>% 
  filter(param == "n") %>% 
  ggplot() + 
  geom_point(aes(iter, true_N), col = "red", size = 2) +
  geom_pointrange(aes(iter, mean, ymin = `2.5%`, ymax = `97.5%`)) +
  ylab("Estimated population") + xlab("Simulation iter.") +
  theme_minimal() +
  coord_flip()


# Check mixing of s, z


temp <- readRDS("intermediate/sim/sim_validation_simplebinomial_uSCR_104.RDS")

plot(as.numeric(unlist(temp$samples[[1]][,"s[100, 1, 1]"])),
     as.numeric(unlist(temp$samples[[1]][,"s[100, 1, 2]"])))

plot(temp$samples[, "z[112]"])
