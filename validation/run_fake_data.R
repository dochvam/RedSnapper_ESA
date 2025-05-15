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

fit_uscr_binom <- function(iter, prefix, M = 500*3, niter = 10000, 
                           nchains = 1, nburnin = 1000, thin = 1, thin2 = 10, 
                           integration_type = "Full") {
  start_time <- Sys.time()
  stopifnot(integration_type %in% c("Full", "Camera_only", "Camera_ROV", 
                                    "Camera_Telemetry", "Camera_only_noCovar"))
  stopifnot(M %% 3 == 0)
  set.seed(14432 + iter * 333) # set seed based on "iter" for reproducibility
  
  source('pipeline_NC/prep_data_NC.R', local = TRUE)
  
  
  #### Randomly generate some counts ####
  y_mtx[] <- rpois(prod(dim(y_mtx)), mean(y_mtx))
  
  n.samples <- sum(y_mtx)
  
  this.j <- this.k <- rep(NA,n.samples)
  idx <- 1
  for(j in 1:J){ #then traps
    for(k in 1:K){ #then occasions
      if(y_mtx[j,k]>0){ #is there at least one sample here?
        for(l in 1:y_mtx[j,k]){ #then samples
          this.j[idx] <- j
          this.k[idx] <- k
          idx <- idx+1
        }
      }
    }
  }
  constants <- list(M = M,
                    J = J,
                    integration_type = integration_type,
                    log_sigma_prior_mean = log_sigma_estimate$mean,
                    log_sigma_prior_sd = log_sigma_estimate$sd,
                    current_dir = camera_locs$current_dir,
                    K1D = rep(K, J),
                    n.samples = n.samples,
                    xlim = xlim,
                    ylim = ylim,
                    datevec = X_datevec,
                    idate = rep(1:3, each = M/3),
                    hm = vps_mtx,
                    ones_mtx = ones_mtx,
                    hm_nrow = nrow(vps_mtx),
                    hm_ncol = ncol(vps_mtx),
                    resoln = res(vps_intensity_ras)[1],
                    nROV = nrow(rov_dat),
                    rb_weights = intersections_df$weight,
                    rov_cell_xvec = intersections_df$x_ind,
                    rov_cell_yvec = intersections_df$y_ind,
                    rbe = rbe, rbs = rbs)
  
  y.true.init <- initialize_ytrue(M,
                                  z_init, 
                                  s_init, 
                                  this.j = this.j,
                                  this.k = this.k,
                                  X_mtx = X_mtx,
                                  current_dir = camera_locs$current_dir,
                                  X_datevec = X_datevec, 
                                  idate = constants$idate,
                                  n.samples = n.samples,
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
                   log_sigma = log_sigma_init,
                   spatial_beta = 0.3)
  
  
  Nimdata <- list(y.true=matrix(NA,nrow=(M),ncol=J),
                  ROV_obs = rov_dat$count,
                  ID = rep(NA, n.samples),
                  z = rep(NA, M),
                  X = as.matrix(X_mtx),
                  capcounts=rep(NA, M))
  
  #### Model code, adapted from Ben Augustine ####
  model_code <- nimbleCode({
    # priors
    lambda.N ~ dunif(0,M*50) #expected abundance
    for (i in 1:3) {
      p0[i] ~ dunif(0,1) #baseline detection probability on logit scale
    }
    
    if (integration_type == "Camera_Telemetry" || integration_type == "Full") {
      log_sigma ~ dnorm(log_sigma_prior_mean, sd = log_sigma_prior_sd) # informative prior around true log sigma
    } else {
      log_sigma ~ dunif(0, 10)
    }
    sigma <- exp(log_sigma)
    
    N ~ dpois(lambda.N) #realized abundance
    
    spatial_beta ~ dnorm(1, sd = 1)
    if (integration_type != "Camera_only_noCovar") {
      phi[1:hm_nrow, 1:hm_ncol] <- exp(spatial_beta * log(hm[1:hm_nrow, 1:hm_ncol])) # hm is log-scale covariate
    } else {
      phi[1:hm_nrow, 1:hm_ncol] <- ones_mtx[1:hm_nrow, 1:hm_ncol]
    }
    
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
                                           current_dir = current_dir[1:J],
                                           idate = idate[i],
                                           p0=p0[1:3], 
                                           z=z[i])
      
      y.true[i,1:J] ~ dBernoulliVector(pd=pd[i, 1:J],
                                       K1D = K1D[1:J], 
                                       z=z[i]) #vectorized obs mod
    }
    #calculate number of inds captured
    capcounts[1:M] <- Getcapcounts(ID=ID[1:n.samples],M=M) #intermediate object
    n <- Getncap(capcounts=capcounts[1:M])
    
    if (integration_type == "Camera_ROV" || integration_type == "Full") {
      for (i in 1:nROV) {
        # rbs and rbe are ROV buffer start/end indexes for ROV i
        pctFishInROVbuffer[i] <- calcPctFishInROVbuffer(phi = phi[1:hm_nrow, 1:hm_ncol], 
                                                        weights = rb_weights[rbs[i]:rbe[i]], 
                                                        rov_cell_xvec = rov_cell_xvec[rbs[i]:rbe[i]],
                                                        rov_cell_yvec = rov_cell_yvec[rbs[i]:rbe[i]],
                                                        n = rbe[i] - rbs[i] + 1)
        ROV_obs[i] ~ dpois(pctFishInROVbuffer[i] * lambda.N * (1/3))
      }
    }
  }) #model
  
  
  #### Build the model ####
  
  parameters <- c('lambda.N', 'p0', 'log_sigma', 'sigma', 'N', 'n', 'spatial_beta')
  
  parameters2 <- c("ID", 's', 'z')
  
  # Build the model, configure the mcmc, and compile
  start.time <- Sys.time()
  Rmodel <- nimbleModel(code=model_code, constants=constants, data=Nimdata,check=FALSE,
                        inits=Niminits)
  
  config.nodes <- c("lambda.N","p0","log_sigma","spatial_beta")
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
                                                    this.j = this.j,
                                                    this.k = this.k,
                                                    n.samples = n.samples,
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
  conf$addSampler(target = c("p0","log_sigma","lambda.N"),
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
               mcmc_time = difftime(mcmc_end_time, mcmc_start_time, units = "mins"),
               total_time = difftime(end_time, start_time, units = "mins"),
               integration_type = integration_type,
               iter = iter,
               prefix = prefix,
               system = "NC_ChickenRock"
  ),
  paste0("intermediate/uSCR_fake_Augustine_Binom", prefix, iter, "_", integration_type, ".RDS"))
}


fit_uscr_binom(iter = 100, prefix = "_TestRandomData_",
               M = 3000, nchains = 3, nburnin = 10000,
               niter = 20000, thin = 2, thin2 = 50,
               integration_type = "Camera_only_noCovar")


result <- readRDS("intermediate/uSCR_fake_Augustine_Binom_TestRandomData_100_Camera_only_noCovar.RDS")
View(result$summary)
