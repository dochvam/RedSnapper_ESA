library(tidyverse)
library(readxl)
library(terra)
library(nimble)

#### Ben Augustine's custom sSampler ####

sSampler <- nimbleFunction(
  # name = 'sampler_RW',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    i <- control$i    
    t <- control$t
    xlim<-control$xlim
    ylim<-control$ylim
    ## control list extraction
    # logScale            <- extractControlElement(control, 'log',                 FALSE)
    # reflective          <- extractControlElement(control, 'reflective',          FALSE)
    adaptive            <- extractControlElement(control, 'adaptive',            TRUE)
    adaptInterval       <- extractControlElement(control, 'adaptInterval',       200)
    adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent', 0.8)
    scale               <- extractControlElement(control, 'scale',               1)
    ## node list generation
    # targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    calcNodes <- model$getDependencies(target)
    # calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
    # isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)   ## should be made faster
    # calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
    # calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
    ## numeric value generation
    scaleOriginal <- scale
    timesRan      <- 0
    timesAccepted <- 0
    timesAdapted  <- 0
    scaleHistory  <- c(0, 0)   ## scaleHistory
    acceptanceHistory  <- c(0, 0)   ## scaleHistory
    if(nimbleOptions('MCMCsaveHistory')) {
      saveMCMChistory <- TRUE
    } else saveMCMChistory <- FALSE
    optimalAR     <- 0.44
    gamma1        <- 0
    ## checks
    # if(length(targetAsScalar) > 1)   stop('cannot use RW sampler on more than one target; try RW_block sampler')
    # if(model$isDiscrete(target))     stop('cannot use RW sampler on discrete-valued target; try slice sampler')
    # if(logScale & reflective)        stop('cannot use reflective RW sampler on a log scale (i.e. with options log=TRUE and reflective=TRUE')
    if(adaptFactorExponent < 0)      stop('cannot use RW sampler with adaptFactorExponent control parameter less than 0')
    if(scale < 0)                    stop('cannot use RW sampler with scale control parameter less than 0')
  },
  run = function() {
    z <- model$z[i]
    if(z==0){#propose from uniform prior
      model$s[i, t, 1:2] <<- c(runif(1, xlim[1], xlim[2]), runif(1, ylim[1], ylim[2]))
      model$calculate(calcNodes)
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    }else{#MH
      s.cand=c(rnorm(1,model$s[i,t,1],scale), rnorm(1,model$s[i,t,2],scale))
      inbox= s.cand[1]< xlim[2] & s.cand[1]> xlim[1] & s.cand[2] < ylim[2] & s.cand[2] > ylim[1]
      if(inbox){
        model_lp_initial <- model$getLogProb(calcNodes)
        model$s[i, t, 1:2] <<- s.cand
        model_lp_proposed <- model$calculate(calcNodes)
        log_MH_ratio <- model_lp_proposed - model_lp_initial
        accept <- decide(log_MH_ratio)
        if(accept) {
          copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        } else {
          copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        }
        if(adaptive){ #we only tune for z=1 proposals
          adaptiveProcedure(accept)
        }
      }
    }
  },
  methods = list(
    adaptiveProcedure = function(jump = logical()) {
      timesRan <<- timesRan + 1
      if(jump)     timesAccepted <<- timesAccepted + 1
      if(timesRan %% adaptInterval == 0) {
        acceptanceRate <- timesAccepted / timesRan
        timesAdapted <<- timesAdapted + 1
        if(saveMCMChistory) {
          setSize(scaleHistory, timesAdapted)                 ## scaleHistory
          scaleHistory[timesAdapted] <<- scale                ## scaleHistory
          setSize(acceptanceHistory, timesAdapted)            ## scaleHistory
          acceptanceHistory[timesAdapted] <<- acceptanceRate  ## scaleHistory
        }
        gamma1 <<- 1/((timesAdapted + 3)^adaptFactorExponent)
        gamma2 <- 10 * gamma1
        adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
        scale <<- scale * adaptFactor
        ## If there are upper and lower bounds, enforce a maximum scale of
        ## 0.5 * (upper-lower).  This is arbitrary but reasonable.
        ## Otherwise, for a poorly-informed posterior,
        ## the scale could grow without bound to try to reduce
        ## acceptance probability.  This creates enormous cost of
        ## reflections.
        # if(reflective) {
        #   lower <- model$getBound(target, 'lower')
        #   upper <- model$getBound(target, 'upper')
        #   if(scale >= 0.5*(upper-lower)) {
        #     scale <<- 0.5*(upper-lower)
        #   }
        # }
        timesRan <<- 0
        timesAccepted <<- 0
      }
    },
    getScaleHistory = function() {       ## scaleHistory
      returnType(double(1))
      if(saveMCMChistory) {
        return(scaleHistory)
      } else {
        print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
        return(numeric(1, 0))
      }
    },          
    getAcceptanceHistory = function() {  ## scaleHistory
      returnType(double(1))
      if(saveMCMChistory) {
        return(acceptanceHistory)
      } else {
        print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
        return(numeric(1, 0))
      }
    },
    ##getScaleHistoryExpanded = function() {                                                 ## scaleHistory
    ##    scaleHistoryExpanded <- numeric(timesAdapted*adaptInterval, init=FALSE)            ## scaleHistory
    ##    for(iTA in 1:timesAdapted)                                                         ## scaleHistory
    ##        for(j in 1:adaptInterval)                                                      ## scaleHistory
    ##            scaleHistoryExpanded[(iTA-1)*adaptInterval+j] <- scaleHistory[iTA]         ## scaleHistory
    ##    returnType(double(1)); return(scaleHistoryExpanded) },                             ## scaleHistory
    reset = function() {
      scale <<- scaleOriginal
      timesRan      <<- 0
      timesAccepted <<- 0
      timesAdapted  <<- 0
      if(saveMCMChistory) {
        scaleHistory  <<- c(0, 0)    ## scaleHistory
        acceptanceHistory  <<- c(0, 0)
      }
      gamma1 <<- 0
    }
  )
)


#### nimbleFunction GetDetectionRate ####

GetDetectionRate <- nimbleFunction(
  run = function(s = double(1), p0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      ans <- p0*exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)

run_one_uSCR_simulation <- function(iter) {
  set.seed(96257 + iter * 41)
  
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
  
  M <- 1000
  psi <- 0.3
  p0 <- 0.9
  log_sigma <- 4
  
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
                                                 sigma = sigma, 
                                                 p0 = p0, 
                                                 z = z[i])
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
    p0 ~ dunif(0, 1)     # Baseline detection
    # log_sigma ~ dnorm(3.435, sd = 1.138) # From telemetry
    log_sigma ~ dnorm(4, sd = 0.1) # Tight prior on true value -- best case
    log(sigma) <- log_sigma
  })
  
  
  mod <- nimbleModel(
    code = my_uscr, 
    constants = list(M = M,
                     ncam = ncam_vec,
                     X = X_array,
                     xlim = xlim,
                     ylim = ylim),
    data = list(
      y = y_mtx
    ),
    inits = list(
      z = rbinom(M, 1, 0.5),
      psi = 0.5,
      p0 = 1,
      log_sigma = 3.5,
      s = array(0, dim = c(M, 3, 2))
    )
  )
  
  cmod <- compileNimble(mod)
  
  mcmcConf <- configureMCMC(cmod)
  mcmcConf$setMonitors(c("p0", "psi", "n", "log_sigma"))
  
  #### Custom samplers
  
  #use block update for p0 and var bc correlated posteriors.
  mcmcConf$removeSampler(c("p0","log_sigma"))
  mcmcConf$addSampler(target = c(paste("p0"),paste("log_sigma")),
                      type = 'AF_slice',control = list(adaptive=TRUE),silent = TRUE)
  
  mcmcConf$removeSampler("s")
  for(i in 1:M){
    for (t in 1:3)
      mcmcConf$addSampler(target = paste("s[",i,",", t, ", 1:2]", sep=""),
                          type = 'sSampler',
                          control=list(i = i, t = t,
                                       xlim = xlim,
                                       ylim = ylim,
                                       scale = 0.25), silent = TRUE)
  }
  
  
  mcmc <- buildMCMC(mcmcConf)
  cmcmc <- compileNimble(mcmc)
  

  res_list <- list()
  
  cmod$p0 <- p0
  cmod$psi  <- psi
  
  cmod$simulate("z")
  cmod$simulate("s")
  true_N <- sum(cmod$z)
  
  cmod$simulate("log_sigma")
  cmod$calculate("sigma")
  cmod$calculate("lam")
  cmod$calculate("lam_sum")
  cmod$simulate("y", includeData = T)
  cmod$setData("y")
  
  samples <- runMCMC(cmcmc, niter = 50000, nburnin = 20000, nchains = 3,
                     thin = 5, samplesAsCodaMCMC = TRUE)
  
  summary <- MCMCvis::MCMCsummary(samples)
  summary$iter <- iter
  summary$param <- rownames(summary)
  summary$true_N <- true_N
  rownames(summary) <- NULL
  
  saveRDS(list(summary = summary,
               samples = samples),
          paste0("intermediate/sim/simple_", iter, ".RDS"))
}
