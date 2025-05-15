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

fit_uscr_binom_oneday <- function(iter, specs_df, prefix, M = 500, niter = 10000, 
                           nchains = 1, nburnin = 1000, thin = 1, thin2 = 10) {
  start_time <- Sys.time()

  integration_type <- specs_df$integration_type[iter]
  target_date <- specs_df$date[iter]
  data_thin_interval <- specs_df$data_thin_interval[iter]
  
  set.seed(999888777 + iter * 333) # set seed based on "iter" for reproducibility
  
  source('pipeline_NC/prep_data_NC.R', local = TRUE)
  
  #### Model code, adapted from Ben Augustine ####
  model_code <- nimbleCode({
    # priors
    lambda.N ~ dunif(0, M*50) #expected abundance
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
        ROV_obs[i] ~ dpois(pctFishInROVbuffer[i] * lambda.N)
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
  
  
  #"sSampler_wCovar", which is a RW block update for the x and y locs with no covariance,
  #and only tuned for when z=1. When z=0, it draws from the distribution on s.
  for(i in 1:(M)){
    conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                    type = 'sSampler_wCovar',
                    control=list(i=i, xlim=xlim, ylim=ylim, scale=50,
                                 resoln = constants$resoln),
                    silent = TRUE)
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
               specs_df_row = specs_df[iter, ],
               system = "NC_ChickenRock"
  ),
  paste0("pipeline_NC/NC_results/oneday/uSCR_OneDay_Binom", prefix, iter, "_", integration_type, ".RDS"))
}



type_vec <- c("Full", "Camera_only", "Camera_ROV", "Camera_Telemetry", "Camera_only_noCovar")
# 
# fit_uscr_binom(iter = 111, prefix = "_LargerBuffer_",
#                M = 3000, nchains = 3, nburnin = 10000,
#                niter = 30000, thin = 2, thin2 = 50,
#                integration_type = "Full")




# Fit models for one day at a time
# Fit models for diff. levels of thinning
# Days 2-3, skip ROV models
# 


specs_df <- expand.grid(
  date = 1:3,
  integration_type = c("Full", "Camera_only", "Camera_ROV", "Camera_Telemetry"),
  data_thin_interval = c(1, 5, 10)
) %>% 
  as.data.frame() %>% 
  filter(date == 1 | integration_type %in% c("Camera_only", "Camera_Telemetry")) %>% 
  mutate(ID = row_number())





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

parLapply(cl, which(specs_df$date == 1 & specs_df$data_thin_interval == 5), fit_uscr_binom_oneday, 
          prefix = "_TestChangingBaseline_", specs_df = specs_df,
          M = 1000, nchains = 3, nburnin = 20000,
          niter = 100000, thin = 2, thin2 = 50)
# parLapply(cl, 1:nrow(specs_df), fit_uscr_binom_oneday, 
#           prefix = "_TestChangingBaseline_", specs_df = specs_df,
#           M = 1000, nchains = 3, nburnin = 20000,
#           niter = 100000, thin = 2, thin2 = 50)




results_files <- list.files("pipeline_NC/NC_results/oneday/",
                            full.names = TRUE)

summary_df <- results_files %>% 
  lapply(function(x) {
    temp <- readRDS(x)
    df <- temp$summary %>% 
      mutate(date = temp$specs_df_row$date,
             thin = temp$specs_df_row$data_thin_interval,
             integration_type = temp$specs_df_row$integration_type,
             ID = temp$specs_df_row$ID
             )
    df$param <- rownames(df)
    rownames(df) <- NULL
    df
  }) %>% 
  bind_rows() 


#### Plot of log(sigma) ####
log_sigma_estimate <- read_csv("pipeline_NC/NC_results/log_sigma_estimate_NC.csv")

summary_df %>% 
  filter(param == "log_sigma") %>% 
  bind_rows(tibble(
    integration_type = "Telem. prior",
    mean = log_sigma_estimate$mean,
    `2.5%` = log_sigma_estimate$mean - 1.96 * log_sigma_estimate$sd,
    `97.5%` = log_sigma_estimate$mean + 1.96 * log_sigma_estimate$sd,
    ID = 0,
    date = 1:3,
    thin = 1
  )) %>% 
  ggplot() + 
  geom_pointrange(aes(ID, mean, ymin = `2.5%`, ymax = `97.5%`,
                      col = integration_type, 
                      shape = factor(paste0("Thin: ", thin, " frame(s)"),
                                     levels = unique(paste0("Thin: ", sort(thin), " frame(s)"))))) +
  facet_wrap(~paste0("Sampling day ", date), ncol = 1) +
  scale_shape_discrete("") +
  coord_flip() +
  theme_bw() + xlab("") +
  theme(axis.ticks = element_blank(), axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  ggtitle("log(sigma)")


#### Plot of N ####
log_sigma_estimate <- read_csv("pipeline_NC/NC_results/log_sigma_estimate_NC.csv")

summary_df %>% 
  filter(param == "N") %>% 
  ggplot() + 
  geom_pointrange(aes(ID, mean, ymin = `2.5%`, ymax = `97.5%`,
                      col = integration_type, 
                      shape = factor(paste0("Thin: ", thin, " frame(s)"),
                                     levels = unique(paste0("Thin: ", sort(thin), " frame(s)"))))) +
  facet_wrap(~paste0("Sampling day ", date), ncol = 1) +
  scale_shape_discrete("") +
  coord_flip() +
  theme_bw() + xlab("") +
  theme(axis.ticks = element_blank(), axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  ggtitle("N")

summary_df %>% 
  filter(param == "spatial_beta") %>% 
  ggplot() + 
  geom_pointrange(aes(ID, mean, ymin = `2.5%`, ymax = `97.5%`,
                      col = integration_type, 
                      shape = factor(paste0("Thin: ", thin, " frame(s)"),
                                     levels = unique(paste0("Thin: ", sort(thin), " frame(s)"))))) +
  facet_wrap(~paste0("Sampling day ", date), ncol = 1) +
  scale_shape_discrete("") +
  coord_flip() +
  theme_bw() + xlab("") +
  theme(axis.ticks = element_blank(), axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  ggtitle("Spatial covar effect")


summary_df %>% 
  filter(param == "spatial_beta") %>% 
  ggplot() + 
  geom_pointrange(aes(ID, mean, ymin = `2.5%`, ymax = `97.5%`,
                      col = integration_type, 
                      shape = factor(paste0("Thin: ", thin, " frame(s)"),
                                     levels = unique(paste0("Thin: ", sort(thin), " frame(s)"))))) +
  facet_wrap(~paste0("Sampling day ", date), ncol = 1) +
  scale_shape_discrete("") +
  coord_flip() +
  theme_bw() + xlab("") +
  theme(axis.ticks = element_blank(), axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  ggtitle("Spatial covar effect")


#### Plot of ESA ####

source("pipeline_NC/ESA_helper.R")
ESA_list <- list() 
for (i in 1:length(results_files)) {
  temp <- readRDS(results_files[i])
  ESA_list[[i]] <- calc_ESA(temp, "Binom") %>% 
    mutate(date = temp$specs_df_row$date,
           thin = temp$specs_df_row$data_thin_interval,
           integration_type = temp$specs_df_row$integration_type,
           ID = temp$specs_df_row$ID)
}


ESA_list %>% 
  bind_rows() %>% 
  ggplot() + 
  geom_pointrange(
    aes(ID, ESA_q50, ymin = ESA_q025, ymax = ESA_q975,
                      col = integration_type, 
                      shape = factor(paste0("Thin: ", thin, " frame(s)"),
                                     levels = unique(paste0("Thin: ", sort(thin), " frame(s)"))))) +
  facet_grid(paste0("Sampling day ", date)~current_dir) +
  scale_shape_discrete("") +
  coord_flip() +
  theme_bw() + xlab("") +
  scale_y_continuous(transform = "log", breaks = 100*2^(2 * 1:5)) +
  theme(axis.ticks = element_blank(), axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  ggtitle("ESA")



