library(tidyverse)
library(readxl)
library(terra)
library(nimble)

if (!dir.exists("intermediate")) dir.create("intermediate")

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

run_one_uSCR_simulation <- function(iter, M = 1000) {
  set.seed(542389 + iter * 41)
  
  #### Load a premade template based on real data ####
  
  if (file.exists("simulations/data_template.RDS")) {
    data_template <- readRDS("simulations/data_template.RDS") 
  } else {
    source("simulations/prep_sim_template_data.R") # <- this file makes the template dat.
  }
  
  #### Set true values for use in simulation ####
  
  psi <- 0.3
  p0 <- 0.9
  log_sigma <- 5
  
  #### Allow for M to change ####
  
  data_template$constants$M <- M
  data_template$inits$z <- rbinom(M, 1, 0.5)
  data_template$inits$s <- array(0, dim = c(M, 3, 2))
  
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
    p0 ~ dgamma(1, 1)     # Baseline detection
    # log_sigma ~ dnorm(3.435, sd = 1.138) # From telemetry
    log_sigma ~ dnorm(5, sd = 0.1) # Tight prior on true value -- best case
    log(sigma) <- log_sigma
  })
  
  
  mod <- nimbleModel(
    code = my_uscr, 
    constants = data_template$constants,
    data = data_template$data,
    inits = data_template$inits
  )
  
  cmod <- compileNimble(mod)
  
  # DON'T add samplers for "s". We want custom samplers and they take
  # forever to remove.
  mcmcConf <- configureMCMC(cmod, nodes = c("psi", "z", "p0", "log_sigma"))
  mcmcConf$setMonitors(c("p0", "psi", "n", "log_sigma"))
  
  #### Custom samplers
  
  #use block update for p0 and var bc correlated posteriors.
  mcmcConf$removeSampler(c("p0","log_sigma"))
  mcmcConf$addSampler(target = c(paste("p0"),paste("log_sigma")),
                      type = 'AF_slice',control = list(adaptive=TRUE),silent = TRUE)
  
  # sampler_targets <- mcmcConf$getSamplers() %>% 
  #   lapply(function(x) x$target) %>% 
  #   unlist()
  
  # mcmcConf$removeSampler(cmod$expandNodeNames("s"))
  for(i in 1:M){
    for (t in 1:3) {
      # mcmcConf$removeSampler(paste("s[",i,",", t, ", 1:2]", sep=""))
      mcmcConf$addSampler(target = paste("s[",i,",", t, ", 1:2]", sep=""),
                          type = 'sSampler',
                          control=list(i = i, t = t,
                                       xlim = data_template$constants$xlim,
                                       ylim = data_template$constants$ylim,
                                       scale = 0.25), silent = TRUE)
    }
  }
  
  
  mcmc <- buildMCMC(mcmcConf)
  cmcmc <- compileNimble(mcmc)
  

  res_list <- list()
  
  cmod$p0 <- p0
  cmod$psi <- psi
  
  cmod$simulate("z")
  cmod$simulate("s")
  true_N <- sum(cmod$z)
  
  cmod$log_sigma <- log_sigma
    
  cmod$calculate("sigma")
  cmod$calculate("lam")
  cmod$calculate("lam_sum")
  cmod$simulate("y", includeData = T)
  cmod$setData("y")
  
  mod$y <- cmod$y
  mod$setData("y")
  
  # 500 iterations w/ M=5000 takes ~15 minutes
  samples <- runMCMC(cmcmc, niter = 50000, nburnin = 20000, nchains = 2,
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
