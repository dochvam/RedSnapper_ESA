library(tidyverse)
library(readxl)
library(nimble)

set.seed(89412)

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
cl <- makeCluster(10)
captured <- clusterEvalQ(cl, source("simulations/SCR_sim_fn.R"))

specs_df <- expand.grid(
  M = c(250, 500), psi = 0.2, log_sigma = log(c(0.1, 0.5, 1, 2.5)),
  grid_edge = c(5, 10, 15, 25), lam0 = c(0.4, 1, 2.5, 6.25), 
  rep = 1:4
) %>% 
  mutate(iter = row_number(), num_detectors = grid_edge^2,
         seed = floor(runif(n()) * 100000))

result <- parLapply(cl = cl, X = 1:nrow(specs_df), 
                    fun = run_one_uSCR_simulation_abstract,
                    prefix = "abstract_big_",
                    specs_df = specs_df)

summary <- list.files("intermediate/sim/", pattern = "abstract_big_", full.names = T) %>% 
  lapply(function(x) readRDS(x)$summary) %>% 
  bind_rows()

summary %>% 
  filter(param == "n") %>% 
  left_join(specs_df, by = "iter") %>% 
  ggplot() + 
  geom_point(aes(as.factor(lam0), true_N, group = iter), col = "red", size = 2,
             position = position_dodge(width = 0.3)) +
  geom_pointrange(aes(as.factor(lam0), mean, group = iter,
                      ymin = `2.5%`, ymax = `97.5%`), position = position_dodge(width = 0.3)) +
  ylab("Population size") + xlab("lambda0") +
  theme_bw() + 
  coord_flip() + 
  theme(axis.ticks = element_blank()) +
  facet_wrap(~factor(paste(num_detectors, "detectors"),
                     levels = paste(sort(unique(specs_df$num_detectors)), "detectors")), scales = "free_y")


summary %>% 
  filter(param == "n") %>% 
  left_join(specs_df, by = "iter") %>% 
  group_by(D = num_detectors, lambda0 = lam0) %>% 
  summarize(abs_error = mean(abs(`50%` - true_N)),
            CI95_width = mean(`97.5%` - `2.5%`)
            # coverage = mean(`97.5%` > true_N & `2.5%` < true_N)
            )
