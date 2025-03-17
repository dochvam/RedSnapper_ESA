source("Ben_Explore/uSCR_fn.R")
library(tidyverse)
library(readxl)
library(terra)

#### Load the camera data ####

camera_dat <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/cameratraps/RS_2023_full_reads_all_three_cameratrap_dates.xlsx") %>% 
  mutate(time = as.numeric(difftime(time, as_datetime("1899-12-31 00:00:00"), units = "sec"))) %>% 
  filter(camera == "A", time >= 60 * 10) # Filter out descent/retrieval frames

camera_counts <- camera_dat %>%
  group_by(Station_ID, date) %>% 
  summarize(count = sum(total))

camera_locs <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/cameratraps/RS_2023_full_reads_all_three_cameratrap_dates.xlsx", sheet = "StationData") %>% 
  filter(`Camera (A or C)` == "A") %>% 
  dplyr::select(Station_ID, Date, Time = Start_Time_GMT,
                Latitude = Start_Latitude, Longitude = Start_Longitude) %>% 
  left_join(camera_counts, by = c("Station_ID", "Date" = "date"))

camera_locs$Longitude <- as.numeric(camera_locs$Longitude)

camera_pts <- vect(camera_locs, geom = c("Longitude", "Latitude"),
                   crs = "+proj=longlat") %>% 
  project("ESRI:102003") %>% 
  as.data.frame(geom = "XY")

X <- camera_pts[, c("x", "y")]
X[, 1] <- X[, 1] - mean(X[, 1])
X[, 2] <- X[, 2] - mean(X[, 2])
xlim <- c(min(X[, 1]) - 500, max(X[, 1]) + 500)
ylim <- c(min(X[, 2]) - 500, max(X[, 2]) + 500)


#### Load hardbottom ####


#### My nimble model ####

my_uscr <- nimbleCode({
  # Loop over all potential individuals
  for (i in 1:M) {
    z[i] ~ dbern(psi)
    
    for (date in 1:3) {
      s[i, 1, date] ~ dunif(xlim[1],xlim[2])
      s[i, 2, date] ~ dunif(ylim[1],ylim[2])
      
      lam[i, 1:ncam, date] <- GetDetectionRate(s = s[i, 1:2, date], 
                                           X = X[1:ncam, 1:2], 
                                           J = ncam, sigma = sigma, 
                                           lam0 = lam0, z = z[i])
    }
  }

  for (k in 1:ncam) {
    lam_sum[k] <- sum(lam[1:M, k, date_vec[k]])
    y[k] ~ dpois(lam_sum[k])
  }
  
  n <- sum(z[1:M])
  
  # Priors
  psi ~ dunif(0,1)     # Inclusion prob. for data augmentation
  lam0 ~ dunif(0,100)   # Baseline detection
  log_var ~ dnorm(7, sd = 7) # Movement param
  sigma <- sqrt(exp(log_var))
})


M <- 2000
mod <- nimbleModel(
  code = my_uscr, 
  constants = list(M = M,
                   ncam = nrow(camera_locs),
                   date_vec = as.numeric(as.factor(camera_locs$Date)),
                   X = as.matrix(X[, c("x", "y")]),
                   xlim = xlim,
                   ylim = ylim
                   ),
  data = list(
    y = camera_locs$count
  ),
  inits = list(
    z = rbinom(M, 1, 0.5),
    psi = 0.5,
    log_var = 10,
    lam0 = 1,
    s = array(0, dim = c(M, 2, 3))
  )
)

cmod <- compileNimble(mod)
mcmcConf <- configureMCMC(cmod)

mcmcConf$setMonitors(c("lam0", "psi", "log_var", "sigma", "n"))


#### Custom samplers

#use block update for lam0 and var bc correlated posteriors.
mcmcConf$removeSampler(c("lam0","log_var"))
mcmcConf$addSampler(target = c(paste("lam0"),paste("log_var")),
                type = 'AF_slice',control = list(adaptive=TRUE),silent = TRUE)

# mcmcConf$removeSampler(paste("s[1:",M,", 1:2, 1:3]", sep=""))
# for(i in 1:M){
#   for (date in 1:3)
#   mcmcConf$addSampler(target = paste("s[",i,", 1:2, ",date,"]", sep=""),
#                   type = 'sSampler', 
#                   control=list(i = i,
#                                xlim = xlim,
#                                ylim = ylim,
#                                scale = 0.25),silent = TRUE)
# }

#### Run the MCMC

mcmc <- buildMCMC(mcmcConf)
cmcmc <- compileNimble(mcmc)

samples <- runMCMC(cmcmc, niter = 10000, nburnin = 5000, thin = 2,
                   nchains = 2, samplesAsCodaMCMC = TRUE)

this_summary <- MCMCvis::MCMCsummary(samples)


plot(samples[, "lam0"])
plot(samples[, "psi"])
plot(samples[, "log_var"])
