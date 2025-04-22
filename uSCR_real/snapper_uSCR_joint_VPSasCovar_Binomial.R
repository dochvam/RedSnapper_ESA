library(nimbleEcology)
library(tidyverse)
library(readxl)
library(terra)
library(parallel)

nimbleOptions(MCMCuseConjugacy = FALSE)
source("uSCR_real/uscr_binom_helper.R")

#### Construct the VPS intensity surface ####
run_one_uscr_chain_binom <- function(iter, prefix, M = 500, niter = 10000, 
                                     nchains = 1, nburnin = 1000, thin = 1,
                                     test_model_only = FALSE) {
  start_time <- Sys.time()
  set.seed(568792 + iter * 333)
  
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
  fish_pts <- vect(fish_positions, geom = c("Longitude", "Latitude"),
                   crs = "+proj=longlat") %>% 
    project(crs(hb_raw))
  
  e <- ext(fish_pts)
  xmin(e) <- floor(xmin(e) / 50) * 50
  xmax(e) <- ceiling(xmax(e) / 50) * 50
  ymin(e) <- floor(ymin(e) / 50) * 50
  ymax(e) <- ceiling(ymax(e) / 50) * 50
  
  template_grid <- rast(e, res = 50)
  crs(template_grid) <- crs(fish_pts)
  terra::values(template_grid) <- 1:ncell(template_grid)
  
  
  cell_counts <- count(extract(template_grid, fish_pts), lyr.1)
  cell_counts <- left_join(data.frame(lyr.1 = 1:ncell(template_grid)), cell_counts)
  cell_counts$n[is.na(cell_counts$n)] <- 0
  cell_counts$prob <- (cell_counts$n+1) / sum(cell_counts$n+1)
  
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
  
  
  #### Load both days of ROV data ####
  
  surface_to_buffer <- vps_intensity_ras 
  source("uSCR_real/process_ROV_buffers.R")
  
  processed_buffers_result <- process_ROV_buffers(buffer_size = 5, surface_to_buffer = surface_to_buffer)
  rbe <- processed_buffers_result$rbe
  rbs <- processed_buffers_result$rbs
  intersections_df <- processed_buffers_result$intersections_df
  rov_dat <- processed_buffers_result$rov_dat
  
  #### Format the camera data ####
  # M <- 1000
  
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
  
  
  
  #### Make model data lists ####
  
  data_template <- list(
    constants = list(M = M,
                     current_dir = current_dir_mtx,
                     ncam = ncam_vec,
                     X = X_array,
                     xmin = grid_bbox[1], 
                     xmax = grid_bbox[2], 
                     ymin = grid_bbox[3], 
                     ymax = grid_bbox[4],
                     resoln = res(vps_intensity_ras)[1], 
                     hm = vps_mtx,
                     hm_nrow = nrow(vps_mtx),
                     hm_ncol = ncol(vps_mtx), 
                     nframes = nframe_mtx,
                     nROV = nrow(rov_dat),
                     # ROV_hb_area = rov_dat$area * rov_dat$structured_pct,
                     rb_weights = intersections_df$weight,
                     rov_cell_xvec = intersections_df$x_ind,
                     rov_cell_yvec = intersections_df$y_ind,
                     rbe = rbe, rbs = rbs
    ),
    data = list(
      y = y_mtx,
      ROV_obs = rov_dat$count
    ),
    inits = list(
      z = matrix(ncol = 3, rbinom(M*3, 1, 0.75)),
      psi = 0.2,
      # lam0 = rep(0.1, 3),
      log_sigma = 4.4,
      sigma = exp(4.4),
      spatial_beta = 0.3,
      p0 = rep(0.1, 3),
      s = initialize_s(n = M, t = 3,
                       xmin = grid_bbox[1],
                       xmax = grid_bbox[2],
                       ymin = grid_bbox[3],
                       ymax = grid_bbox[4],
                       resoln = 50, 
                       habitatMask = vps_mtx,
                       spatial_beta = 0.4)
    )
  )
  
  
  #### nimble model code ####
  
  my_uscr <- nimbleCode({
    
    for (i in 1:nROV) {
      # rbs and rbe are ROV buffer start/end indexes for ROV i
      pctFishInROVbuffer[i] <- calcPctFishInROVbuffer(phi = phi[1:hm_nrow, 1:hm_ncol], 
                                                      weights = rb_weights[rbs[i]:rbe[i]], 
                                                      rov_cell_xvec = rov_cell_xvec[rbs[i]:rbe[i]],
                                                      rov_cell_yvec = rov_cell_yvec[rbs[i]:rbe[i]],
                                                      n = rbe[i] - rbs[i] + 1)
      ROV_obs[i] ~ dpois(pctFishInROVbuffer[i] * psi * M)
    }
    
    phi[1:hm_nrow, 1:hm_ncol] <- exp(spatial_beta * log(hm[1:hm_nrow, 1:hm_ncol])) # hm is log-scale covariate
    
    # Loop over all potential individuals
    for (i in 1:M) {
      # Latent state representing the inclusion prob.
      
      for (t in 1:3) {
        z[i, t] ~ dbern(psi)
        
        # Distribution of centroids
        s[i, t, 1:2] ~ dHabDistr_asCovar(
          xmax = xmax,
          xmin = xmin,
          ymax = ymax,
          ymin = ymin,
          resoln = resoln,
          phi = phi[1:hm_nrow, 1:hm_ncol]
        )
        
        # Calculate distances
        detprob[i, 1:ncam[t], t] <- GetDetectionRate(s = s[i, t, 1:2], 
                                                    X = X[1:ncam[t], t, 1:2], 
                                                    J = ncam[t], 
                                                    sigma = sigma, 
                                                    lam0 = p0[1:3],
                                                    current_dir = current_dir[1:ncam[t], t],
                                                    z = z[i, t])
      }
    }
    
    for (t in 1:3) {
      for (j in 1:ncam[t]) {
        y[j, t] ~ dPoisBinom_wReps(detprob[1:M, j, t], reps = nframes[j, t]) # per frame
      }
      
      n[t] <- sum(z[1:M, t])
    }
    
    
    # Priors
    psi ~ dbeta(0.01, 1) # Scale prior for inclusion prob.
    for (i in 1:3) {
      p0[i] ~ dunif(0, 1)  # Detection rate per frame at centroid, depends on current dir.
    }
    spatial_beta ~ dnorm(1, sd = 1)
    
    log_sigma ~ dnorm(3.435, sd = 1.138) # From telemetry
    log(sigma) <- log_sigma
    
    # sigma ~ dunif(1, 1000) # <- uninformative prior
    # log_sigma <- log(sigma)
  })
  
  
  mod <- nimbleModel(
    code = my_uscr, 
    constants = data_template$constants,
    data = data_template$data,
    inits = data_template$inits,
    calculate = F
  )
  
  cmod <- compileNimble(mod)
  if (test_model_only) return(cmod$calculate())
  
  
  mcmcConf <- configureMCMC(cmod, nodes = c("s", "spatial_beta"))
  mcmcConf$addSampler(target = c("p0", "log_sigma", "psi"), type = "AF_slice")

  # combos <- expand.grid(1:3, 1:2)
  nodeControl <- list(mean = 0, scale = 1)
  for (i in 1:M) {
    for (t in 1:3) {
      nodeControl$targetNode_c <- paste0("s[", i, ", ", t, ", 1:2]") # supply concatenated version of variable to deal with model[[coefNode]]
      nodeControl$targetNode <- paste0("s[", i, ",", t, ",", 1:2, "]")
      mcmcConf$addSampler(type = my_sampler_MV_RJ_indicator,
                          target = paste0("z[", i, ",", t, "]"),
                          control = nodeControl)

      this_node <- paste0("s[", i, ",", t, ", 1:2]")
      ## Add sampler for the coefficient variable (when is in the model)
      currentConf <- mcmcConf$getSamplers(this_node)
      mcmcConf$removeSamplers(this_node)
      mcmcConf$addSampler(type = my_sampler_RJ_toggled,
                          target = this_node,
                          control = list(samplerType = currentConf[[1]],
                                         fixedValue = c(0,0)))
    }
  }
  
  # 
  # # We want custom samplers for everything except z, s
  # mcmcConf <- configureMCMC(cmod, nodes = c("z", "spatial_beta"))
  # 
  # #### Custom samplers
  # #use block update for lam0, psi, and var bc correlated posteriors.
  # mcmcConf$addSampler(target = c("lam0", "psi", "log_sigma"),
  #                     # mcmcConf$addSampler(target = c("lam0", "psi", "sigma"),
  #                     type = 'AF_slice', control = list(adaptive=TRUE), silent = TRUE)
  
  mcmcConf$setMonitors(c("p0", "psi", "n", "log_sigma", "sigma", "spatial_beta", "s", "z"))
  
  
  mcmc <- buildMCMC(mcmcConf)
  cmcmc <- compileNimble(mcmc)
  
  mcmc_start_time <- Sys.time()
  samples <- runMCMC(cmcmc, niter = niter, nburnin = nburnin, nchains = nchains,
                     thin = thin, samplesAsCodaMCMC = TRUE, 
                     inits = data_template$inits)
  mcmc_end_time <- Sys.time()
  
  # samples <- runMCMC(cmcmc, niter = 2000, nburnin = 0, nchains = 1,
  #                    thin = 1, samplesAsCodaMCMC = TRUE)
  
  summary <- MCMCvis::MCMCsummary(samples)
  summary$param <- rownames(summary)
  
  rownames(summary) <- NULL
  end_time <- Sys.time()
  
  saveRDS(list(summary = summary,
               samples = samples,
               mcmc_time = difftime(mcmc_end_time, mcmc_start_time, units = "mins"),
               total_time = difftime(end_time, start_time, units = "mins")
               ),
          # paste0("uSCR_real/joint_masked_VPSasCovar_Pois_uninformativePrior.RDS"))
          paste0("uSCR_real/", prefix, iter, "_Binom.RDS"))
}




cl <- makeCluster(5)
capture <- clusterEvalQ(cl, {
  library(nimbleEcology)
  library(tidyverse)
  library(readxl)
  library(terra)
  library(parallel)
  
  nimbleOptions(MCMCuseConjugacy = FALSE)
  source("uSCR_real/uscr_binom_helper.R")
})

parLapply(cl, 1:5, run_one_uscr_chain_binom, 
          prefix = "uSCR_real_RJMCMC_", niter = 100, M = 500, 
          nburnin = 0, thin = 1, nchains = 1, test_model_only = FALSE
          )

stopCluster(cl)
rm(cl)


summary_all <- lapply(list.files("uSCR_real/", pattern = "uSCR_real_RJ", full.names = TRUE),
       function(x) {
         readRDS(x)$summary %>% 
           mutate(fn = x)
       }) %>% 
  bind_rows()
times <- lapply(list.files("uSCR_real/", pattern = "uSCR_real_RJ", full.names = TRUE),
       function(x) {
         readRDS(x)$mcmc_time
       }) %>% 
  unlist()
