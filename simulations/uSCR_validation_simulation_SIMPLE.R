library(tidyverse)
library(readxl)
library(nimble)

##########################################################
#' Do we get anything useful out of our SCR model at the
#' scale we're working at?
#' 
#' In this simulation:
#'  > Constant sigma w/ informed prior based on 6-hr interval
#'  > Ignoring time-staggered deployment starts
#'  > Real camera locations
##########################################################
source("simulations/SCR_sim_fn.R")

#### Parallel execution ####
library(parallel)
cl <- makeCluster(8)
captured <- clusterEvalQ(cl, source("simulations/SCR_sim_fn.R"))

result <- parLapply(cl = cl, X = 8 + 1:8, fun = run_one_uSCR_simulation,
                    M = 5000)

summary <- list.files("intermediate/sim/", pattern = "simple_", full.names = T) %>% 
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


