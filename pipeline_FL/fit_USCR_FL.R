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

fit_uscr_binom_FL <- function(iter, prefix, M = 500*2, niter = 10000, 
                           nchains = 1, nburnin = 1000, thin = 1, thin2 = 10, 
                           spatial_type_vec, integration_type_vec, sigma_type_vec,
                           M_type_vec) {
  start_time <- Sys.time()
  
  spatial_type <- spatial_type_vec[iter]
  integration_type <- integration_type_vec[iter]
  sigma_type <- sigma_type_vec[iter]
  
  stopifnot(spatial_type     %in% c("HB_map", "VPS_map", "SpConstant"))
  stopifnot(sigma_type       %in% c("Prior_Mean", "Prior_Variability", "Trunc_prior", "Weak_prior", "Constant"))
  stopifnot(integration_type %in% c("Full_Offset", "Full_NoOffset", "No_ROV"))
  stopifnot(M_type_vec       %in% c("MBaseline", "MAugment"))
  
  if (M_type_vec[iter] == "MAugment") M <- M * 1.5
  
  stopifnot(M %% 2 == 0)
  set.seed(14432 + iter * 222) # set seed based on "iter" for reproducibility
  
  source('pipeline_FL/prep_data_FL.R', local = TRUE)
  
  #### Model code, adapted from Ben Augustine ####
  model_code <- nimbleCode({
    # priors
    lambda.N ~ dgamma(1e-6, 1e-6) #expected abundance
    for (i in 1:3) {
      p0[i] ~ dunif(0,1) #baseline detection probability on logit scale
    }
    
    if (sigma_type == "Prior_Mean" || sigma_type == "Prior_Variability") {
      log_sigma_long ~ dnorm(log_sigma_prior_mean, sd = log_sigma_prior_sd) # informative prior around true log sigma
    } else if (sigma_type == "Trunc_prior") {
      log_sigma_long ~ T(dnorm(log_sigma_prior_mean, sd = log_sigma_prior_sd), 
                         log_sigma_prior_mean - 3*log_sigma_prior_sd, 
                         log_sigma_prior_mean + 3*log_sigma_prior_sd)
    } else if (sigma_type == "Weak_prior") {
      log_sigma_long ~ dunif(0, 10)
    }else if (sigma_type == "Constant"){
      log_sigma_long <- log_sigma_prior_mean
    }
    
    sigma <- exp(log_sigma_long)
    
    N ~ dpois(lambda.N) #realized abundance

    if (spatial_type == "SpConstant") {
      phi[1:hm_nrow, 1:hm_ncol] <- exp(0 * log(hm[1:hm_nrow, 1:hm_ncol])) 
    } else {
      phi[1:hm_nrow, 1:hm_ncol] <- exp(spatial_beta * log(hm[1:hm_nrow, 1:hm_ncol])) 
    }
      spatial_beta ~ dnorm(1, sd = 1)
    
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
    
    if (integration_type == "Full_NoOffset" || integration_type == "Full_Offset") {
      for (i in 1:nROV) {
        # rbs and rbe are ROV buffer start/end indexes for ROV i
        pctFishInROVbuffer[i] <- calcPctFishInROVbuffer(phi = phi[1:hm_nrow, 1:hm_ncol], 
                                                        weights = rb_weights[rbs[i]:rbe[i]], 
                                                        rov_cell_xvec = rov_cell_xvec[rbs[i]:rbe[i]],
                                                        rov_cell_yvec = rov_cell_yvec[rbs[i]:rbe[i]],
                                                        n = rbe[i] - rbs[i] + 1)
        ROV_obs[i] ~ dpois(pctFishInROVbuffer[i] * lambda.N * (1/3) * ROV_offset)
      }
    }
    if (integration_type == "Full_Offset") {
      ROV_offset ~ dunif(0, 10)
    } else {
      ROV_offset <- 1
    }
  }) #model
  
  
  #### Build the model ####
  
  if (sigma_type == "Constant") {
    parameters <- c('lambda.N', 'p0','log_sigma_long', 'sigma', 'N', 'n', 'spatial_beta')
    config.nodes <- c("lambda.N", "p0","spatial_beta")
  } else {
    parameters <- c('lambda.N', 'p0', 'log_sigma_long', 'sigma', 'N', 'n', 'spatial_beta')
    config.nodes <- c("lambda.N","p0","log_sigma_long","spatial_beta")
  }
  
  if (integration_type == "Full_Offset") {
    parameters <- c(parameters, "ROV_offset")
    config.nodes <- c(config.nodes, "ROV_offset")
  } 
  
  parameters2 <- c("ID", 's', 'z')
  
  # Build the model, configure the mcmc, and compile
  start.time <- Sys.time()
  Rmodel <- nimbleModel(code=model_code, constants=constants, data=Nimdata,check=FALSE,
                        inits=Niminits)
  
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
  
  
  for(i in 1:(M)){
    conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                    type = 'sSampler_wCovar',
                    control=list(i=i, xlim=xlim, ylim=ylim, scale=50,
                                 resoln = constants$resoln),
                    silent = TRUE)
    #scale parameter here is just the starting scale. It will be tuned.
  }
  
  #use block update for  correlated posteriors. Can use "tries" to control how many times per iteration
  if (sigma_type == "Constant") {
    conf$addSampler(target = c("p0","lambda.N"),
                    type = 'RW_block',control = list(adaptive=TRUE,tries=1),silent = TRUE)
  } else {
    conf$addSampler(target = c("p0","log_sigma_long","lambda.N"),
                    type = 'RW_block',control = list(adaptive=TRUE,tries=1),silent = TRUE)
  }
  
  
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
               spatial_type = spatial_type,
               sigma_type = sigma_type,
               M_type = M_type_vec[iter],
               iter = iter,
               prefix = prefix,
               system = "FL_TurtleMount"
  ),
  paste0("pipeline_FL/FL_results/uSCR_real_Augustine_Binom", prefix, iter, "_", 
         integration_type, "_", sigma_type, "_", spatial_type, "_", M_type_vec[iter], ".RDS"))
}

cl <- makeCluster(8)

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
  
  nimbleOptions(determinePredictiveNodesInModel = FALSE)
  
})


combos <- expand.grid(
  spatial_type  = c("VPS_map"),
  integration_type = c("Full_Offset", "No_ROV"),
  sigma_type = c("Prior_Mean", "Prior_Variability", "Weak_prior"),
  M_type = c("MBaseline", "MAugment")
) %>% 
  as.data.frame()

nchains <- 4

parLapply(cl, 1:(nrow(combos)*nchains), fit_uscr_binom_FL, prefix = "_FinalChoices_",
          spatial_type_vec = rep(combos$spatial_type, nchains),
          integration_type_vec = rep(combos$integration_type, nchains),
          sigma_type_vec = rep(combos$sigma_type, nchains),
          M_type_vec = rep(combos$M_type, nchains),
          M = 15000, nchains = 1, nburnin = 40000,
          niter = 100000, thin = 1, thin2 = 20)

# parLapply(cl, 1:6, fit_uscr_binom_FL, prefix = "_TEST_",
#           spatial_type_vec = rep(combos$spatial_type, 1),
#           integration_type_vec = rep(combos$integration_type, 1),
#           sigma_type_vec = rep(combos$sigma_type, 1),
#           M = 10000, nchains = 1, nburnin = 0,
#           niter = 1000, thin = 1, thin2 = 1)


# fit_uscr_binom(iter = 100, prefix = "_Final_",
#                M = 3000, nchains = 3, nburnin = 100000,
#                niter = 200000, thin = 2, thin2 = 50,
#                integration_type = "Full")

stopCluster(cl)
