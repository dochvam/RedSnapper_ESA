### In this script---use Ben Augustine's code to implement the binomial USCR.
###    Use a simulation framework, but try to add our real locations and
###    also add the VPS surface.
library(tidyverse)
library(readxl)
library(terra)
library(parallel)
library(nimble)
library(coda)
library(MCMCvis)

# Source file: Ben Augustine's helper fns
source("uSCR_binom_Augustine/augustine_helper.R")

# Source file: Ben Goldstein's helper fns
source("uSCR_binom_Augustine/other_helper.R")

nimbleOptions(determinePredictiveNodesInModel = FALSE)

binomial_sim_per_Augustine <- function(iter, prefix, M = 500*3, niter = 10000, 
                                       nchains = 1, nburnin = 1000, thin = 1, thin2 = 10) {
  start_time <- Sys.time()
  
  stopifnot(M %% 3 == 0)
  
  start_time <- Sys.time()
  set.seed(568792 + iter * 333) # set seed based on "iter" for reproducibility
  
  #### Process the real data ####
  # First, we process all the real data. This allows us to define simulation data
  # dimensions, use the real VPS data (as spatial covar.) and camera locations,
  # etc. We ultimately only simulate the observed counts at each camera and ROV
  
  hb_raw <- rast("Data/Chicken_Rock_Map/ChickenRock_Classification.tif")
  
  # Get all the VPS fixes of living fish
  VPS_folder <- "Data/SnapperMvmtAbundanceStudy/VPS_Data/VPS-ChickenRock-01-Results-20240202/results/animal/"
  VPS_files <- list.files(VPS_folder)
  fish_files <- VPS_files[grepl("^12", VPS_files)]
  fish_fate <- read_xlsx("Data/SnapperMvmtAbundanceStudy/VPS_Data/Fate_assignments_from_discard_paper_vps 2023.xlsx")
  fish_to_drop <- fish_fate$`Tag number`[fish_fate$`subsequent assignments based on velocity and depth data` == "Discard mortality"]
  dead_fish_files <- paste0(fish_to_drop, ".csv")
  fish_files <- fish_files[!fish_files %in% dead_fish_files]
  fish_positions <- lapply(file.path(VPS_folder, fish_files),
                           read_csv) %>% 
    bind_rows()
  
  # Convert to pts
  fish_pts <- vect(fish_positions, geom = c("Longitude", "Latitude"),
                   crs = "+proj=longlat") %>% 
    project(crs(hb_raw))
  
  # Define a spatial extent encompassing all VPS locs, snapped to a 50 m grid
  e <- ext(fish_pts)
  xmin(e) <- floor(xmin(e) / 50) * 50
  xmax(e) <- ceiling(xmax(e) / 50) * 50
  ymin(e) <- floor(ymin(e) / 50) * 50
  ymax(e) <- ceiling(ymax(e) / 50) * 50
  
  template_grid <- rast(e, res = 50)
  crs(template_grid) <- crs(fish_pts)
  terra::values(template_grid) <- 1:ncell(template_grid)
  
  # Count the number of VPS fixes in each grid cell, then get a total proportion
  cell_counts <- count(extract(template_grid, fish_pts), lyr.1)
  cell_counts <- left_join(data.frame(lyr.1 = 1:ncell(template_grid)), cell_counts)
  cell_counts$n[is.na(cell_counts$n)] <- 0
  # Add 1 obs. to every cell so we don't make it impossible for a fish to be somewhere
  cell_counts$prob <- (cell_counts$n+1) / sum(cell_counts$n+1)
  
  # Formatting stuff
  vps_intensity_ras <- template_grid
  terra::values(vps_intensity_ras) <- cell_counts$prob
  plot(vps_intensity_ras)
  
  vps_mtx <- t(as.matrix(vps_intensity_ras, wide = T))
  vps_mtx <- vps_mtx[, ncol(vps_mtx):1]
  
  grid_bbox <- ext(vps_intensity_ras)
  grid_bbox <- as.numeric(c(grid_bbox[1], grid_bbox[2], grid_bbox[3], grid_bbox[4]))
  
  x_offset <- mean(grid_bbox[1:2])
  y_offset <- mean(grid_bbox[3:4])
  
  grid_bbox[1:2] <- grid_bbox[1:2] - x_offset
  grid_bbox[3:4] <- grid_bbox[3:4] - y_offset
  
  resoln <- res(vps_intensity_ras)[1]
  
  
  #### Load both days of ROV data ####
  
  surface_to_buffer <- vps_intensity_ras
  source("uSCR_real/process_ROV_buffers.R")
  
  processed_buffers_result <- process_ROV_buffers(buffer_size = 4.9, surface_to_buffer = surface_to_buffer)
  rbe <- processed_buffers_result$rbe
  rbs <- processed_buffers_result$rbs
  intersections_df <- processed_buffers_result$intersections_df
  rov_dat <- processed_buffers_result$rov_dat
  
  #### Load and format the camera data ####
  
  # Get the real camera locations and deployment times for use in the simulation
  
  camera_dat <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/cameratraps/RS_2023_full_reads_all_three_cameratrap_dates.xlsx") %>% 
    mutate(time = as.numeric(difftime(time, as_datetime("1899-12-31 00:00:00"), units = "sec"))) %>% 
    filter(camera == "A", time > 60 * 10, time <= 60 * 30) # Filter out descent/retrieval frames
  
  camera_counts <- camera_dat %>%
    group_by(Station_ID, date) %>% 
    summarize(count = sum(total), nframes = n())
  
  camera_locations_corrected <- read_csv("intermediate/corrected_camera_stations.csv") %>% 
    select(Station_ID, Date, Longitude, Latitude)
  
  camera_locs <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/cameratraps/RS_2023_full_reads_all_three_cameratrap_dates.xlsx", sheet = "StationData") %>% 
    filter(`Camera (A or C)` == "A") %>%
    dplyr::select(Station_ID, Date, Time = Start_Time_GMT, 
                  current_dir = `Current Direction`) %>% 
    left_join(camera_counts, by = c("Station_ID", "Date" = "date")) %>% 
    left_join(camera_locations_corrected, by = c("Station_ID", "Date"))
  
  camera_locs$current_dir <- ifelse(
    grepl("towards", tolower(camera_locs$current_dir)), 1, 
    ifelse(grepl("away", tolower(camera_locs$current_dir)), 2, 3)
  )
  
  
  camera_locs$Longitude <- as.numeric(camera_locs$Longitude)
  
  camera_pts <- vect(camera_locs, geom = c("Longitude", "Latitude"),
                     crs = "+proj=longlat") %>% 
    project(crs(hb_raw)) %>% 
    as.data.frame(geom = "XY")
  
  X <- camera_pts[, c("x", "y")]
  X[, 1] <- X[, 1] - x_offset
  X[, 2] <- X[, 2] - y_offset
  
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
  
  current_dir_mtx <- y_mtx 
  for (i in 1:3) {
    ncam_vec[i] <- sum(camera_locs$Date == Dates[i])
    current_dir_mtx[1:ncam_vec[i], i] <- camera_locs$current_dir[camera_locs$Date == Dates[i]]
  }
  
  
  # Make a 3D array of camera coords  
  X_array <- array(dim = c(mtx_nrow, 3, 2))
  for (i in 1:3) {
    X_array[1:ncam_vec[i], i, 1] <- X$x[camera_locs$Date == Dates[i]]
    X_array[1:ncam_vec[i], i, 2] <- X$y[camera_locs$Date == Dates[i]]
  }
  # Convert to a 2d matrix
  X_mtx <- rbind(X_array[,1,], X_array[,2,], X_array[1:14,3,])
  X_datevec <- c(rep(1, 18), rep(2, 18), rep(3, 14))
  stopifnot(nrow(X_mtx) == length(X_datevec))
  
  #### Simulate observed counts in an R function ####
  
  # Simulate camera data #
  xlim <- grid_bbox[1:2]
  ylim <- grid_bbox[3:4]
  J <- nrow(X_mtx)
  
  K <- 20
  true_spatial_beta <- 0.3
  true_psi <- 0.5
  true_p0 <- 0.8
  true_sigma <- exp(4)
  
  sim.out <- sim_counts_wHabCovar(
    M = M, psi = true_psi, p0 = true_p0, sigma = true_sigma, 
    X_mtx = X_mtx, X_datevec = X_datevec,
    J = J, xlim = xlim, ylim = ylim, K = K,
    resoln = resoln,
    habitatMask = vps_mtx,
    spatial_beta = true_spatial_beta
  )
  
  # this.j and this.k describe the samplers and occasions on which each detection occurred
  
  # Simulate ROV data #
  ROV_counts <- sim_ROV(nROV = nrow(rov_dat),
                        rb_weights = intersections_df$weight,
                        rov_cell_xvec = intersections_df$x_ind,
                        rov_cell_yvec = intersections_df$y_ind,
                        rbe = rbe, rbs = rbs, 
                        spatial_beta = true_spatial_beta, 
                        habitatMask = vps_mtx)
  
  #### Make NIMBLE model input lists ####
  
  constants <- list(M = M,
                    J = J,
                    K1D = rep(K, J),
                    n.samples = sim.out$n.samples,
                    xlim = xlim,
                    ylim = ylim,
                    datevec = X_datevec,
                    idate = rep(1:3, each = M/3),
                    hm = vps_mtx,
                    hm_nrow = nrow(vps_mtx),
                    hm_ncol = ncol(vps_mtx),
                    resoln = res(vps_intensity_ras)[1],
                    nROV = nrow(rov_dat),
                    # ROV_hb_area = rov_dat$area * rov_dat$structured_pct,
                    rb_weights = intersections_df$weight,
                    rov_cell_xvec = intersections_df$x_ind,
                    rov_cell_yvec = intersections_df$y_ind,
                    rbe = rbe, rbs = rbs)
  
  
  z_init <- rbinom(M, 1, 0.5)
  s_init <- initialize_s(
    M = M,
    xlim = xlim,
    ylim = ylim,
    resoln = resoln,
    habitatMask = vps_mtx,
    spatial_beta = 0.5
  )
  sigma_init <- exp(4.5)
  p0_init <- 0.9
  y.true.init <- initialize_ytrue(M,
                                  z_init, 
                                  s_init, 
                                  this.j = sim.out$this.j,
                                  this.k = sim.out$this.k,
                                  X_datevec = X_datevec, 
                                  idate = constants$idate,
                                  n.samples = sim.out$n.samples,
                                  sigma_init = sigma_init, p0_init = p0_init, K = K)
  
  Niminits <- list(z = z_init,
                   N = sum(z_init), #must initialize N to be the sum of z init
                   lambda.N=sum(z_init), #initializing lambda.N to be consistent with N.init
                   s = s_init,
                   ID = y.true.init$ID,
                   capcounts = rowSums(y.true.init$ytrue2D),
                   y.true = y.true.init$ytrue2D,
                   p0 = p0_init, 
                   sigma = sigma_init,
                   spatial_beta = 0.3)
  
  
  Nimdata <- list(y.true=matrix(NA,nrow=(M),ncol=J),
                  ROV_obs = ROV_counts,
                  ID = rep(NA, sim.out$n.samples),
                  z = rep(NA, M),
                  X = as.matrix(X_mtx),
                  capcounts=rep(NA, M))
  
  #### Model code, adapted from Ben Augustine ####
  model_code <- nimbleCode({
    # priors
    lambda.N ~ dunif(0,M*2) #expected abundance
    logit(p0) ~ dlogis(0,1) #baseline detection probability on logit scale
    # sigma ~ dunif(0,10) #detection spatial scale parameter, uninformative
    sigma ~ dnorm(exp(3.5), 5) #informative
    N ~ dpois(lambda.N) #realized abundance
    spatial_beta ~ dnorm(1, sd = 1)
    
    phi[1:hm_nrow, 1:hm_ncol] <- exp(spatial_beta * log(hm[1:hm_nrow, 1:hm_ncol])) # hm is log-scale covariate
    
    for(i in 1:M) {
      
      s[i, 1:2] ~ dHabDistr_asCovar(
        xmax = xlim[2],
        xmin = xlim[1],
        ymax = ylim[2],
        ymin = ylim[1],
        resoln = resoln,
        phi = phi[1:hm_nrow, 1:hm_ncol]
      )

      pd[i,1:J] <- GetDetectionProb_wDates(s = s[i,1:2], 
                                           X = X[1:J, 1:2], 
                                           J=J, 
                                           sigma=sigma, 
                                           datevec = datevec[1:J],
                                           idate = idate[i],
                                           p0=p0, 
                                           z=z[i])
      
      y.true[i,1:J] ~ dBernoulliVector(pd=pd[i, 1:J],
                                       K1D = K1D[1:J], 
                                       z=z[i]) #vectorized obs mod
    }
    #calculate number of inds captured
    capcounts[1:M] <- Getcapcounts(ID=ID[1:n.samples],M=M) #intermediate object
    n <- Getncap(capcounts=capcounts[1:M])
    
    for (i in 1:nROV) {
      # rbs and rbe are ROV buffer start/end indexes for ROV i
      pctFishInROVbuffer[i] <- calcPctFishInROVbuffer(phi = phi[1:hm_nrow, 1:hm_ncol], 
                                                      weights = rb_weights[rbs[i]:rbe[i]], 
                                                      rov_cell_xvec = rov_cell_xvec[rbs[i]:rbe[i]],
                                                      rov_cell_yvec = rov_cell_yvec[rbs[i]:rbe[i]],
                                                      n = rbe[i] - rbs[i] + 1)
      ROV_obs[i] ~ dpois(pctFishInROVbuffer[i] * lambda.N * (1/3))
    }
  }) #model
  
  
  #### Build the model ####
  
  parameters <- c('lambda.N', 'p0', 'sigma', 'N', 'n', 'spatial_beta')

  parameters2 <- c("ID", 's', 'z')

  # Build the model, configure the mcmc, and compile
  start.time <- Sys.time()
  Rmodel <- nimbleModel(code=model_code, constants=constants, data=Nimdata,check=FALSE,
                        inits=Niminits)
  
  config.nodes <- c("lambda.N","logit_p0","sigma","spatial_beta")
  # config.nodes <- c()
  conf <- configureMCMC(Rmodel,monitors=parameters, thin=thin, 
                        monitors2=parameters2, thin2=thin2, nodes=config.nodes,
                        useConjugacy = FALSE) 
  
  #conf$printSamplers() #shows the samplers used for each parameter and latent variable
  ###Two *required* sampler replacements
  ##Here, we remove the default sampler for y.true
  #and replace it with the custom "IDSampler".
  # conf$removeSampler("y.true")
  #how many times to propose updates for each sample ID per iteration. No idea what is optimal in specific scenarios.
  IDups <- 2
  conf$addSampler(target = paste0("y.true[1:",M,",1:",J,"]"),
                  type = 'IDSampler',control = list(M = M, J = J, K=K,
                                                    this.j = sim.out$this.j,
                                                    this.k = sim.out$this.k,
                                                    n.samples = sim.out$n.samples,
                                                    IDups = IDups),
                  silent = TRUE)
  
  z.ups <- round(M*0.25) # how many N/z proposals per iteration? Not sure what is optimal, setting to 50% of M here.
  #nodes used for update
  y.nodes <- Rmodel$expandNodeNames(paste("y.true[1:",M,",1:",J,"]"))
  pd.nodes <- Rmodel$expandNodeNames(paste("pd[1:",M,",1:",J,"]"))
  N.node <- Rmodel$expandNodeNames(paste("N"))
  z.nodes <- Rmodel$expandNodeNames(paste("z[1:",M,"]"))
  calcNodes <- c(N.node,pd.nodes,y.nodes)
  conf$addSampler(target = c("N"),
                  type = 'zSampler',control = list(z.ups=z.ups,M=M,J=J,
                                                   y.nodes=y.nodes,pd.nodes=pd.nodes,
                                                   N.node=N.node,z.nodes=z.nodes,
                                                   calcNodes=calcNodes),silent = TRUE)
  
  
  #"sSampler", which is a RW block update for the x and y locs with no covariance,
  #and only tuned for when z=1. When z=0, it draws from the prior, assumed to be uniform. 
  # conf$removeSampler(paste("s[1:",M,", 1:2]", sep=""))
  # BRG note: atm we let the sampler propose draws from a uniform prior even though
  #  the actual distribution on s is not uniform. I don't think this is an issue--
  #  just inefficient--but I wouldn't mind a santy check on this
  for(i in 1:(M)){
    conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                    type = 'sSampler',control=list(i=i,xlim=xlim,ylim=ylim,scale=50),silent = TRUE)
    #scale parameter here is just the starting scale. It will be tuned.
  }
  
  #use block update for  correlated posteriors. Can use "tries" to control how many times per iteration
  conf$addSampler(target = c("logit_p0","sigma","lambda.N"),
                  type = 'RW_block',control = list(adaptive=TRUE,tries=1),silent = TRUE)
  
  
  # Build and compile
  Rmcmc <- buildMCMC(conf)
  # runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
  Cmodel <- compileNimble(Rmodel)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  
  
  mcmc_start_time <- Sys.time()
  # Run the MCMC, conventionally
  mcmc_samples <- runMCMC(
    Cmcmc, niter = niter, nburnin = nburnin, thin = thin, thin2 = thin2, nchains = nchains,
    # Cmcmc, niter = 1000, nburnin = 200, thin = 1, thin2 = 50, nchains = 1,
    samplesAsCodaMCMC = TRUE
  )
  mcmc_end_time <- Sys.time()
  
  # Get a summary df of the main parameters
  summary <- MCMCvis::MCMCsummary(mcmc_samples$samples)
  
  end_time <- Sys.time()
  
  # Save the output
  saveRDS(list(summary = summary,
               samples = mcmc_samples$samples,
               samples2 = mcmc_samples$samples2,
               truth_df = data.frame(
                 p0 = true_p0,
                 spatial_beta = true_spatial_beta,
                 sigma = true_sigma,
                 N = sim.out$trueN,
                 N3 = sim.out$trueN3
               ),
               mcmc_time = difftime(mcmc_end_time, mcmc_start_time, units = "mins"),
               total_time = difftime(end_time, start_time, units = "mins")
  ),
  paste0("intermediate/sim/uSCR_Augustine_Binom", prefix, iter, ".RDS"))
}


# Run 5 chains in parallel
cl <- makeCluster(5)

# Load all packages and helper fns within each process environment
capture <- clusterEvalQ(cl, {
  library(tidyverse)
  library(readxl)
  library(terra)
  library(parallel)
  library(nimble)
  library(coda)
  library(MCMCvis)
  
  # Source file: Ben Augustine's helper fns
  source("uSCR_binom_Augustine/augustine_helper.R")
  
  # Source file: Ben Goldstein's helper fns
  source("uSCR_binom_Augustine/other_helper.R")
})

# Run one simulation chains in parallel
parLapply(cl, 10 + 1:5, run_one_uscr_chain_binom, 
          prefix = "_sim_", niter = 50000, M = 1500, 
          nburnin = 5000, thin = 1, thin2 = 25, nchains = 2)

stopCluster(cl)
rm(cl)
