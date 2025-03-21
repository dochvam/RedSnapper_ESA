
#### nimbleFunction GetDetectionRate ####
GetDetectionRate <- nimbleFunction(
  run = function(s = double(1), lam0=double(0), sigma=double(0), 
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


#### Format the camera data ####
M <- 1000

# Get the real camera locations and deployment times for use in the simulation

camera_dat <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/cameratraps/RS_2023_full_reads_all_three_cameratrap_dates.xlsx") %>% 
  mutate(time = as.numeric(difftime(time, as_datetime("1899-12-31 00:00:00"), units = "sec"))) %>% 
  filter(camera == "A", time > 60 * 10, time <= 60 * 30) # Filter out descent/retrieval frames

camera_counts <- camera_dat %>%
  group_by(Station_ID, date) %>% 
  summarize(count = sum(total), nframes = n())

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
xlim <- c(min(X[, 1]) - 100, max(X[, 1]) + 100)
ylim <- c(min(X[, 2]) - 100, max(X[, 2]) + 100)

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

# Make a matrix of frames per camera
nframe_mtx <- y_mtx
for (i in 1:3) {
  ncam_vec[i] <- sum(camera_locs$Date == Dates[i])
  nframe_mtx[1:ncam_vec[i], i] <- camera_locs$nframes[camera_locs$Date == Dates[i]]
}

# Make a 3D array of camera coords  
X_array <- array(dim = c(mtx_nrow, 3, 2))
for (i in 1:3) {
  X_array[1:ncam_vec[i], i, 1] <- X$x[camera_locs$Date == Dates[i]]
  X_array[1:ncam_vec[i], i, 2] <- X$y[camera_locs$Date == Dates[i]]
}



data_template <- list(
  constants = list(M = M,
                   ncam = ncam_vec,
                   X = X_array,
                   xlim = xlim,
                   ylim = ylim,
                   nframes = nframe_mtx),
  data = list(
    y = y_mtx
  ),
  inits = list(
    z = matrix(ncol = 3, rbinom(M*3, 1, 0.1)),
    psi = 0.1,
    lam0 = 0.1,
    log_sigma = 3,
    s = array(0, dim = c(M, 3, 2))
  )
)

#### nimble model code ####

my_uscr <- nimbleCode({
  # Loop over all potential individuals
  for (i in 1:M) {
    # Latent state representing the inclusion prob.
    
    for (t in 1:3) {
      z[i, t] ~ dbern(psi)
      
      # Distribution of centroids
      s[i, t, 1] ~ dunif(xlim[1],xlim[2])
      s[i, t, 2] ~ dunif(ylim[1],ylim[2])
      
      # Calculate distances
      lambda[i, 1:ncam[t], t] <- GetDetectionRate(s = s[i, t, 1:2], 
                                                  X = X[1:ncam[t], t, 1:2], 
                                                  J = ncam[t], 
                                                  sigma = sigma, 
                                                  lam0 = lam0, 
                                                  z = z[i, t])
    }
  }
  
  for (t in 1:3) {
    for (j in 1:ncam[t]) {
      lambda_sum[j, t] <- sum(lambda[1:M, j, t] * nframes[j, t])
      y[j, t] ~ dpois(lambda_sum[j, t])
    }
    
    n[t] <- sum(z[1:M, t])
  }
  
  
  # Priors
  psi ~ dbeta(0.01, 1) # Scale prior for inclusion prob.
  lam0  ~ dgamma(0.1, 1)  # Detection rate per frame at centroid
  
  log_sigma ~ dnorm(3.435, sd = 1.138) # From telemetry
  log(sigma) <- log_sigma
})


mod <- nimbleModel(
  code = my_uscr, 
  constants = data_template$constants,
  data = data_template$data,
  inits = data_template$inits
)

cmod <- compileNimble(mod)

# We want custom samplers for everything except z, s
mcmcConf <- configureMCMC(cmod, nodes = c("z", "s"))

#### Custom samplers
#use block update for lam0 and var bc correlated posteriors.
mcmcConf$addSampler(target = c("lam0", "psi", "log_sigma"),
                    type = 'AF_slice', control = list(adaptive=TRUE), silent = TRUE)

mcmcConf$setMonitors(c("lam0", "psi", "n", "log_sigma", "sigma"))


mcmc <- buildMCMC(mcmcConf)
cmcmc <- compileNimble(mcmc)


samples <- runMCMC(cmcmc, niter = 30000, nburnin = 5000, nchains = 8,
                   thin = 2, samplesAsCodaMCMC = TRUE, 
                   inits = data_template$inits)
# samples <- runMCMC(cmcmc, niter = 2000, nburnin = 0, nchains = 1,
#                    thin = 1, samplesAsCodaMCMC = TRUE)

summary <- MCMCvis::MCMCsummary(samples)
summary$param <- rownames(summary)

rownames(summary) <- NULL

saveRDS(list(summary = summary,
             samples = samples),
        paste0("uSCR_real/simple_posterior.RDS"))
