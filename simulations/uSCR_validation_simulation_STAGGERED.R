library(tidyverse)
library(readxl)

#### nimbleFunction GetDetectionRate ####

GetDetectionRate <- nimbleFunction(
  run = function(s = double(1), lam0=double(0), sigma=double(1), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      ans <- lam0*exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)


#### Load the camera data ####

# Get the real camera locations and deployment times for use in the simulation

camera_dat <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/cameratraps/RS_2023_full_reads_all_three_cameratrap_dates.xlsx") %>% 
  mutate(time = as.numeric(difftime(time, as_datetime("1899-12-31 00:00:00"), units = "sec"))) %>% 
  filter(camera == "A", time >= 60 * 10) # Filter out descent/retrieval frames

camera_counts <- camera_dat %>%
  group_by(Station_ID, date) %>% 
  summarize(count = max(total))

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

# For each camera, get "time since first deployment began"
camera_locs <- camera_locs %>% 
  group_by(Date) %>% 
  mutate(mins_since_start = as.numeric(difftime(Time, min(Time), units = "mins")) + 20)

Dates <- unique(camera_locs$Date)

#### Data reformatting for model ####

# Make a matrix of observed counts
mtx_nrow <- max(table(camera_locs$Date))
ncam_vec <- numeric(3)
y_mtx <- matrix(ncol = 3, nrow = mtx_nrow)
for (i in 1:3) {
  ncam_vec[i] <- sum(camera_locs$Date == Dates[i])
  y_mtx[1:ncam_vec[i], i] <- camera_locs$count[camera_locs$Date == Dates[i]]
}


# Make a matrix of sigma_mu
sigma_mu_mtx <- matrix(ncol = 3, nrow = mtx_nrow)
for (i in 1:3) {
  sigma_mu_mtx[1:ncam_vec[i], i] <- 
    0.431 * log(camera_locs$mins_since_start[camera_locs$Date == Dates[i]]) + 1.002
}

X_array <- array(dim = c(mtx_nrow, 3, 2))
for (i in 1:3) {
  X_array[1:ncam_vec[i], i, 1] <- X$x[camera_locs$Date == Dates[i]]
  X_array[1:ncam_vec[i], i, 2] <- X$y[camera_locs$Date == Dates[i]]
}

#### Set true values for use in simulation ####

M <- 500
psi <- 0.2
lam0 <- 0.9

#### My nimble model ####

my_uscr <- nimbleCode({
  
  # Loop over all potential individuals
  for (i in 1:M) {
    # Latent state representing the inclusion prob.
    z[i] ~ dbern(psi)
    
    for (t in 1:3) {
      # Distribution of centroids
      s[i, t, 1] ~ dunif(xlim[1],xlim[2])
      s[i, t, 2] ~ dunif(ylim[1],ylim[2])
      
      # Calculate distances
      lam[i, 1:ncam[t], t] <- GetDetectionRate(s = s[i, t, 1:2], 
                                                     X = X[1:ncam[t], t, 1:2], 
                                                     J = ncam[t], 
                                                     sigma = sigma[i, 1:ncam[t], t], 
                                                     lam0 = lam0, 
                                                     z = z[i])
      
      # Distribution of random movement parameters sigma
      for (j in 1:ncam[t]) {
        # Save some time by calculating sigma_mu outside of model
        log_sigma[i, j, t] ~ dnorm(sigma_mu[j, t], sd = 1.138)
        log(sigma[i, j, t]) <- log_sigma[i, j, t]
      }
    }
  }
  
  for (t in 1:3) {
    for (j in 1:ncam[t]) {
      lam_sum[j, t] <- sum(lam[1:M, j, t])
      y[j, t] ~ dpois(lam_sum[j, t])
    }
  }
  
  n <- sum(z[1:M])
  
  # Priors
  psi  ~ dunif(0, 1)   # Inclusion prob. for data augmentation
  lam0 ~ dunif(0, 1)   # Baseline detection
})


mod <- nimbleModel(
  code = my_uscr, 
  constants = list(M = M,
                   ncam = ncam_vec,
                   X = X_array,
                   sigma_mu = sigma_mu_mtx,
                   xlim = xlim,
                   ylim = ylim
  ),
  data = list(
    y = y_mtx
  ),
  inits = list(
    z = rbinom(M, 1, 0.5),
    psi = 0.5,
    lam0 = 1,
    log_sigma = array(1, dim = c(M, max(ncam_vec), 3)),
    s = array(0, dim = c(M, 3, 2))
  )
)

cmod <- compileNimble(mod)

mcmcConf <- configureMCMC(cmod)
mcmcConf$setMonitors(c("lam0", "psi", "n"))

#### Custom samplers

#use block update for lam0 and var bc correlated posteriors.
# mcmcConf$removeSampler(c("lam0","log_var"))
# mcmcConf$addSampler(target = c(paste("lam0"),paste("log_var")),
#                     type = 'AF_slice',control = list(adaptive=TRUE),silent = TRUE)

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


mcmc <- buildMCMC(mcmcConf)
cmcmc <- compileNimble(mcmc)

set.seed(96257)

res_list <- list()
for (sim in 1:nsim) {
  cmod$lam0 <- lam0
  cmod$psi  <- psi
  
  cmod$simulate("z")
  cmod$simulate("s")
  
  cmod$simulate("log_sigma")
  cmod$calculate("sigma")
  cmod$calculate("lam")
  cmod$calculate("lam_sum")
  cmod$simulate("y", includeData = T)
  cmod$setData("y")
  
  samples <- runMCMC(cmcmc, niter = 1000, nburnin = 500, nchains = 2,
                     samplesAsCodaMCMC = TRUE)
  
  summary <- MCMCvis::MCMCsummary(samples)
  summary$iter <- sim
  summary$param <- rownames(summary)
  rownames(summary) <- NULL
  res_list[[sim]] <- summary
}


samples <- runMCMC(cmcmc, niter = 10000, nburnin = 5000, thin = 2,
                   nchains = 2, samplesAsCodaMCMC = TRUE)

this_summary <- MCMCvis::MCMCsummary(samples)


plot(samples[, "lam0"])
plot(samples[, "psi"])
plot(samples[, "log_var"])
