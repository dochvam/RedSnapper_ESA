library(nimble)

#------------------------------------------------------------------
# Function for calculation detection rate
#------------------------------------------------------------------
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

dPoissonVector <- nimbleFunction(
  run = function(x = double(1), lam = double(1), z = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      logProb <- sum(dpois(x, lambda = lam, log = TRUE))
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rPoissonVector <- nimbleFunction(
  run = function(n = integer(0),lam = double(1),z = double(0)) {
    returnType(double(1))
    J <- nimDim(lam)[1]
    out <- numeric(J,value=0)
    return(out)
  }
)

Getcapcounts <- nimbleFunction(
  run = function(y.true=double(2)){
    returnType(double(1))
    M <- nimDim(y.true)[1]
    J <- nimDim(y.true)[2]
    capcounts <- numeric(M, value = 0)
    for(i in 1:M){
      capcounts[i] <- sum(y.true[i,1:J])
    }
    return(capcounts)
  }
)
Getncap <- nimbleFunction(
  run = function(capcounts=double(1),ID=double(1)){ #don't need ID, but nimble requires is it used in a function 
    returnType(double(0))
    M <- nimDim(capcounts)[1]
    nstate <- numeric(M, value = 0)
    for(i in 1:M){
      if(capcounts[i]>0){
        nstate[i] <- 1
      }
    }
    n.cap <- sum(nstate)
    return(n.cap)
  }
)

## sampler to update y[1:M,1:J]
IDSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    this.j <- control$this.j
    trapup <- control$trapup
    j.indicator <- control$j.indicator
    M <- control$M
    J <- control$J
    cluster.ups <- control$cluster.ups
    local.eval <- control$local.eval
    swap.rad.multiplier <- control$swap.rad.multiplier
    calcNodes <- model$getDependencies(c("y.true","z"))
  },
  
  run = function() {
    z <- model$z
    for(j in 1:length(trapup)){ #only iterate through traps with samples
      lam.curr <- model$lam[,trapup[j]] #individual by trap expected counts
      lam.curr[z==0] <- 0 #can zero out z==0 here, will already be zeroed out using supplied nimble functions
      lam.use <- lam.curr
      sum.lam.use=sum(lam.use)
      if(sum.lam.use>0){ #abort if 0. Should never happen?
        #full conditional for identity update at this j
        fullcond <- lam.use/sum.lam.use
        idx <- which(this.j==trapup[j]) #which samples were recorded at this trap?
        for(l in 1:length(idx)){ #update the individual identities of these samples one by one
          #update derived parameter ID
          ID.curr <- model$ID[idx[l]]
          ID.cand <- rcat(1,fullcond)
          model$ID[idx[l]] <<- ID.cand
          #update y
          model$y.true[ID.curr,trapup[j]] <<- model$y.true[ID.curr,trapup[j]]-1 #subtract out old ID obs
          model$y.true[ID.cand,trapup[j]] <<- model$y.true[ID.cand,trapup[j]]+1 #add in new ID obs
        }
      }
    }
    
    #Now, we do a joint z-ID update
    if(cluster.ups>0){ #skip if cluster.ups=0
      y.true <- model$y.true
      lam <- model$lam
      K1D <- model$K1D
      z <- model$z
      #precalculate ll.y
      ll.y <- matrix(0,M,J)
      for(i in 1:M){
        if(z[i]==1){
          ll.y[i,] = dpois(y.true[i,],K1D*lam[i,],log=TRUE)
        }
      }
      ID.curr2 <- model$ID #can't reuse object with same name but different dimensions, adding "2" to make nimble happy
      swap.rad <- swap.rad.multiplier*model$sigma[1] #radius for traps to update sample IDs around a focal individual
      for(up in 1:cluster.ups){ #how many times to do this update per iteration?
        # select z==1 to turn off
        z.cand <- z
        ID.cand2 <- ID.curr2
        y.cand <- y.true
        ll.y.cand <- ll.y
        lam.cand <- lam
        z.on <- which(z==1)
        n.z.on <- length(z.on)
        if(n.z.on>1){ #Cannot turn off anyone unless there are at least 2 guys on. samples must belong to someone!
          if(n.z.on>1){
            pick <- rcat(1,rep(1/n.z.on,n.z.on))
            focal.i <- z.on[pick]
          }else{
            focal.i <- z.on[1]
          }
          z.cand[focal.i] <- 0
          lam.cand[focal.i,] <- 0
          p.select.z.for <- 1/n.z.on
          if(local.eval==TRUE){# find local traps with samples
            dists <- sqrt((model$s[focal.i,1]-model$X[,1])^2+(model$s[focal.i,2]-model$X[,2])^2)
            focal.traps <- which(dists<swap.rad&j.indicator) #j.indicator removes traps with 0 samples
          }else{
            focal.traps <- which(j.indicator) #j.indicator removes traps with 0 samples
          }
          total.log.j.probs.for <- 0
          total.log.j.probs.back <- 0
          n.focal.traps <- length(focal.traps)
          abort <- FALSE #abort if any propprobs so small we get underflow. Would be rejected if there were no underflow.
          #almost never happens...
          if(n.focal.traps>0){
            # repropose all samples at focal.traps
            for(j in 1:n.focal.traps){
              these.samps <- which(this.j==focal.traps[j])
              n.these.samps <- length(these.samps)
              propprobs.for <- lam.cand[,focal.traps[j]]*z.cand
              propprobs.back <- lam[,focal.traps[j]]*z
              sum.propprobs.for <- sum(propprobs.for)
              if(sum.propprobs.for==0){
                abort <- TRUE
              }
              propprobs.for <- propprobs.for/sum.propprobs.for
              propprobs.back <- propprobs.back/sum(propprobs.back)
              y.cand[,focal.traps[j]] <- 0
              for(l in 1:n.these.samps){
                pick <- rcat(1,prob=propprobs.for)
                ID.cand2[these.samps[l]] <- pick
                y.cand[ID.cand2[these.samps[l]],focal.traps[j]] <- y.cand[ID.cand2[these.samps[l]],focal.traps[j]]+1
              }
              total.log.j.probs.for <- total.log.j.probs.for+dmulti(y.cand[,focal.traps[j]],n.these.samps,prob=propprobs.for,log=TRUE)
              total.log.j.probs.back <- total.log.j.probs.back+dmulti(y.true[,focal.traps[j]],n.these.samps,prob=propprobs.back,log=TRUE)
              
              #update ll.y.cand - only focal traps with samples here
              for(i in 1:M){
                if(z.cand[i]==1){
                  ll.y.cand[i,focal.traps[j]] <- dpois(y.cand[i,focal.traps[j]],
                                                       K1D[focal.traps[j]]*lam.cand[i,focal.traps[j]],log=TRUE)
                }
              }
            }
            #update ll.y.cand for focal.i
            ll.y.cand[focal.i,] <- 0
          }else{#if we only turn off a z and no local samples to reallocate
            ll.y.cand[focal.i,] <- 0
          }
          if(!abort){#if propprobs didn't have underflow
            ll.z.curr <- dbinom(z[focal.i],1,model$psi[1],log=TRUE)
            ll.z.cand <- dbinom(z.cand[focal.i],1,model$psi[1],log=TRUE)
            
            z.off <- which(z.cand==0)
            p.select.z.back <- 1/length(z.off)
            
            logforprob <- log(p.select.z.for)+total.log.j.probs.for
            logbackprob <- log(p.select.z.back)+total.log.j.probs.back
            
            if(n.focal.traps>0){#y likelihood of everyone at focal traps and all traps for focal individual
              ll.total.curr <- sum(ll.y)+ll.z.curr #just summing full y likelihood for ease
              ll.total.cand <- sum(ll.y.cand)+ll.z.cand
            }else{#y likelihood for focal ind only, all traps
              ll.total.curr <- sum(ll.y[focal.i,])+ll.z.curr
              ll.total.cand <- sum(ll.y.cand[focal.i,])+ll.z.cand
            }
            log_MH_ratio <- (ll.total.cand+logbackprob)-(ll.total.curr+logforprob)
            accept <- decide(log_MH_ratio)
            if(accept){
              if(n.focal.traps>0){
                y.true[,focal.traps] <- y.cand[,focal.traps]
                ll.y[,focal.traps] <- ll.y.cand[,focal.traps]
                ID.curr2 <- ID.cand2
              }
              lam[focal.i,] <- lam.cand[focal.i,]
              ll.y[focal.i,] <- ll.y.cand[focal.i,]
              z[focal.i] <- z.cand[focal.i]
            }
          }
        }
        
        #select z==0 to turn on
        z.cand <- z
        ID.cand2 <- ID.curr2
        y.cand <- y.true
        ll.y.cand <- ll.y
        lam.cand <- lam
        z.off <- which(z==0)
        n.z.off <- length(z.off)
        if(n.z.off>0){
          if(n.z.off>1){
            pick <- rcat(1,rep(1/n.z.off,n.z.off))
            focal.i <- z.off[pick]
          }else{
            focal.i <- z.off[1]
          }
          z.cand[focal.i] <- 1
          
          p.select.z.for <- 1/length(z.off)
          
          dists <- sqrt((model$s[focal.i,1]-model$X[,1])^2+(model$s[focal.i,2]-model$X[,2])^2)
          if(local.eval==TRUE){# find local traps with samples
            focal.traps <- which(dists<swap.rad&j.indicator) #j.indicator removes traps with 0 samples
          }else{
            focal.traps <- which(j.indicator) #j.indicator removes traps with 0 samples
          }
          lam.cand[focal.i,] <- model$lam0[1]*exp(-dists^2/(2*model$sigma[1]^2))
          total.log.j.probs.for <- 0
          total.log.j.probs.back <- 0
          n.focal.traps <- length(focal.traps)
          if(n.focal.traps>0){
            #repropose all samples at focal.traps
            for(j in 1:n.focal.traps){
              these.samps <- which(this.j==focal.traps[j])
              n.these.samps <- length(these.samps)
              propprobs.for <- lam.cand[,focal.traps[j]]*z.cand
              propprobs.back <- lam[,focal.traps[j]]*z
              propprobs.for <- propprobs.for/sum(propprobs.for)
              propprobs.back <- propprobs.back/sum(propprobs.back)
              y.cand[,focal.traps[j]] <- 0
              for(l in 1:n.these.samps){
                pick <- rcat(1,prob=propprobs.for)
                ID.cand2[these.samps[l]] <- pick
                y.cand[ID.cand2[these.samps[l]],focal.traps[j]] <- y.cand[ID.cand2[these.samps[l]],focal.traps[j]]+1
              }
              total.log.j.probs.for <- total.log.j.probs.for+dmulti(y.cand[,focal.traps[j]],n.these.samps,prob=propprobs.for,log=TRUE)
              total.log.j.probs.back <- total.log.j.probs.back+dmulti(y.true[,focal.traps[j]],n.these.samps,prob=propprobs.back,log=TRUE)
              #update ll.y.cand - only focal traps with samples here
              for(i in 1:M){
                if(z.cand[i]==1){
                  ll.y.cand[i,focal.traps[j]] <- dpois(y.cand[i,focal.traps[j]],
                                                       K1D[focal.traps[j]]*lam.cand[i,focal.traps[j]],log=TRUE)
                }
              }
            }
            ll.y.cand[focal.i,] <- dpois(y.cand[focal.i,],K1D*lam.cand[focal.i,],log=TRUE)
          }else{#if we only turn on a z and no local samples to reallocate
            ll.y.cand[focal.i,] <- dpois(y.cand[focal.i,],K1D*lam.cand[focal.i,],log=TRUE)
          }
          ll.z.curr <- dbinom(z[focal.i],1,model$psi[1],log=TRUE)
          ll.z.cand <- dbinom(z.cand[focal.i],1,model$psi[1],log=TRUE)
          
          z.on <- which(z.cand==1)
          p.select.z.back <- 1/length(z.on)
          
          logforprob <- log(p.select.z.for)+total.log.j.probs.for
          logbackprob <- log(p.select.z.back)+total.log.j.probs.back
          
          if(n.focal.traps>0){#y likelihood of everyone at focal traps and all traps for focal individual
            ll.total.curr <- sum(ll.y)+ll.z.curr #just summing full likelihood for ease
            ll.total.cand <- sum(ll.y.cand)+ll.z.cand
          }else{#y likelihood for focal ind only, all traps
            ll.total.curr <- sum(ll.y[focal.i,])+ll.z.curr
            ll.total.cand <- sum(ll.y.cand[focal.i,])+ll.z.cand
          }
          
          log_MH_ratio <- (ll.total.cand+logbackprob)-(ll.total.curr+logforprob)
          accept <- decide(log_MH_ratio)
          if(accept){
            if(n.focal.traps>0){
              y.true[,focal.traps] <- y.cand[,focal.traps]
              ID.curr2 <- ID.cand2
              ll.y[,focal.traps] <- ll.y.cand[,focal.traps]
            }
            ll.y[focal.i,] <- ll.y.cand[focal.i,]
            z[focal.i] <- z.cand[focal.i]
            lam[focal.i,] <- lam.cand[focal.i,]
          }
        }
      }
      
      #update model$stuff
      model$y.true <<- y.true
      model$ID <<- ID.curr2
      model$z <<- z #lam will be updated with calculate below
    }
    #update lp
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)

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
    z <- model$z[i, t]
    if(z==0){#propose from uniform prior
      model$s[i, 1:2] <<- c(runif(1, xlim[1], xlim[2]), runif(1, ylim[1], ylim[2]))
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

e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}


simUSCR <-
  function(N=120,lam0=NA,p0=NA,sigma=0.50,theta=NA,lambda=NA,
           K=10,X=X,buff=3,obstype="poisson"){
    # simulate a population of activity centers
    s <- cbind(runif(N, min(X[,1])-buff, max(X[,1])+buff), 
               runif(N, min(X[,2])-buff, max(X[,2])+buff))
    D <- e2dist(s,X)
    J <- nrow(X)
    
    # Capture individuals
    y <- array(0,dim=c(N,J,K))
    if(obstype=="bernoulli"){
      if(is.na(p0))stop("must provide p0 for bernoulli obstype")
      pd <- p0*exp(-D*D/(2*sigma*sigma))
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y[i,j,k] <- rbinom(1,1,pd[i,j])
          }
        }
      }
    }else if(obstype=="poisson"){
      if(is.na(lam0))stop("must provide lam0 for poisson obstype")
      lamd <- lam0*exp(-D*D/(2*sigma*sigma))
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y[i,j,k] <- rpois(1,lamd[i,j])
          }
        }
      } 
    }else if(obstype=="negbin"){
      if(is.na(lam0))stop("must provide lam0 for negbin obstype")
      if(is.na(theta))stop("Must provide theta for negbin obstype")
      lamd <- lam0*exp(-D*D/(2*sigma*sigma))
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y[i,j,k] <- rnbinom(1,mu=lamd[i,j],size=theta)
          }
        }
      } 
    }else if(obstype=="hurdleZTPois"){
      if(is.na(p0))stop("must provide p0 for hurdleZTpois obstype")
      if(is.na(lambda))stop("must provide lambda for hurdleZTpois obstype")
      library(VGAM)
      pd<- p0*exp(-D*D/(2*sigma*sigma))
      y.det=array(0,dim=c(N,J,K))
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y.det[i,j,k] <- rbinom(1,1,pd[i,j])
            if(y.det[i,j,k]==1){
              y[i,j,k] <- rzapois(1,lambda,pobs0=0)
            }
          }
        }
      }
    }else{
      stop("obstype not recognized")
    }
    
    #discard uncaptured inds and aggregate true IDcovs for all samples, keeping track of where they came from with A matrix (used to put observed data back together)
    caught <- which(apply(y,c(1),sum)>0)
    y.true <- y
    y <- y[caught,,]
    if(K==1){
      y <- array(y,dim=c(dim(y),1))
    }
    n <- length(caught)
    n.samples <- sum(y)
    
    ID <- this.j <- this.k <- rep(NA,n.samples)
    idx <- 1
    for(i in 1:length(caught)){ #loop through inds (uncaptured already removed)
      for(j in 1:J){ #then traps
        for(k in 1:K){ #then occasions
          if(y[i,j,k]>0){ #is there at least one sample here?
            for(l in 1:y[i,j,k]){ #then samples
              ID[idx] <- i
              this.j[idx] <- j
              this.k[idx] <- k
              idx <- idx+1
            }
          }
        }
      }
    }
    
    #reconstruct y to make sure this algorithm works
    y.check <- y*0
    for(l in 1:n.samples){
      y.check[ID[l],this.j[l],this.k[l]] <- y.check[ID[l],this.j[l],this.k[l]]+1
    }
    if(!all(y==y.check))stop("Error rebuilding data. Report bug.")
    
    out <- list(this.j=this.j,this.k=this.k,#observed data
                ID=ID,y=y,y.true=y.true,X=X,K=K,buff=buff,obstype=obstype,s=s,n=nrow(y))
    return(out)
  }


#### nimbleCode ####
NimModel <- nimbleCode({
  #--------------------------------------------------------------
  # priors
  #--------------------------------------------------------------
  psi ~ dunif(0,1)
  lam0 ~ dunif(0,10)
  sigma ~ dunif(0,100)
  #--------------------------------------------------------------
  for(i in 1:M) {
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, lam0=lam0, z=z[i])
    y.true[i,1:J] ~ dPoissonVector(lam=lam[i,1:J]*K1D[1:J],z=z[i]) #vectorized obs mod
  }
  # capcounts[1:M] <- Getcapcounts(y.true=y.true[1:M,1:J])
  # n <- Getncap(capcounts=capcounts[1:M],ID=ID[1:n.samples])
  N <- sum(z[1:M])
})#model

e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.data.USCR <- function(data=NA,M=NA,inits=inits,obstype="poisson"){
  library(abind)
  this.j <- data$this.j
  this.k <- data$this.k
  X <- as.matrix(data$X)
  J <- nrow(X)
  K <- data$K
  n.samples <- length(this.j)
  
  #state space extent
  buff <- data$buff
  xlim <- c(min(X[,1]),max(X[,1]))+c(-buff, buff)
  ylim <- c(min(X[,2]),max(X[,2]))+c(-buff, buff)
  
  # psi <- inits$psi
  sigma <- inits$sigma
  
  #assign random activity centers
  s <- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2]))
  z <- rep(1,M)
  D <- e2dist(s, X)
  
  if(obstype%in%c("poisson","negbin")){
    lam0 <- inits$lam0
    lamd <- lam0*exp(-D*D/(2*sigma*sigma))
    #Build y.true
    y.true <- array(0,dim=c(M,J,K))
    ID <- rep(NA,n.samples)
    for(l in 1:n.samples){
      propdist <- z*lamd[,this.j[l]]
      propdist <- propdist/sum(propdist)
      ID[l] <- sample(1:M,1,replace=FALSE,prob=propdist)
      y.true[ID[l],this.j[l],this.k[l]] <- y.true[ID[l],this.j[l],this.k[l]]+1
    }
  }else if(obstype=="bernoulli"){
    p0 <- inits$p0
    pd <- p0*exp(-D*D/(2*sigma*sigma))
    #Build y.true
    y.true <- array(0,dim=c(M,J,K))
    ID <- rep(NA,n.samples)
    for(l in 1:n.samples){
      propdist <- z*pd[,this.j[l]]
      propdist <- propdist*(1-y.true[,this.j[l],this.k[l]]) #zero out ID's that already have a sample at this j-k
      propdist <- propdist/sum(propdist)
      ID[l] <- sample(1:M,1,replace=FALSE,prob=propdist)
      y.true[ID[l],this.j[l],this.k[l]] <- y.true[ID[l],this.j[l],this.k[l]]+1
    }
  }else if(obstype=="hurdleZTPois"){
    p0 <- inits$p0
    pd <- p0*exp(-D*D/(2*sigma*sigma))
    y.true <- array(0,dim=c(M,J,K))
    ID <- rep(NA,n.samples)
    for(l in 1:n.samples){
      propdist <- z*pd[,this.j[l]]
      propdist <- propdist/sum(propdist)
      ID[l] <- sample(1:M,1,replace=FALSE,prob=propdist)
      y.true[ID[l],this.j[l],this.k[l]] <- y.true[ID[l],this.j[l],this.k[l]]+1
    }
  }else if(obstype=="ramsey"){
    p0 <- inits$p0
    pd <- p0*exp(-D*D/(2*sigma*sigma))
    #Build y.true
    y.true <- array(0,dim=c(M,J,K))
    jk.idx <- which(data$y.jk==1,arr.ind=TRUE)
    n.j.k <- nrow(jk.idx)
    #just allocate 1 individual detection where there was at least 1
    for(l in 1:n.j.k){
      propdist <- z*pd[,jk.idx[l,1]]
      this.i <- which(propdist==max(propdist))
      y.true[this.i,jk.idx[l,1],jk.idx[l,2]] <- 1
    }
  }
  
  y.true2D <- apply(y.true,c(1,2),sum)
  z <- 1*(rowSums(y.true2D)>0)
  
  #Optimize s after assigning samples
  idx <- which(rowSums(y.true2D)>0) #switch for those actually caught
  for(i in idx){
    trps <- matrix(X[y.true2D[i,]>0,1:2],ncol=2,byrow=FALSE)
    if(nrow(trps)>1){
      s[i,] <- c(mean(trps[,1]),mean(trps[,2]))
    }else{
      s[i,] <- trps
    }
  }
  
  sigma <- inits$sigma
  D <- e2dist(s, X)
  if(obstype=="bernoulli"|obstype=="ramsey"){
    p0 <- inits$p0
    pd <- p0*exp(-D*D/(2*sigma*sigma))
    ll.y <- dbinom(y.true2D,K,pd*z,log=TRUE)
  }else if(obstype=="poisson"){
    lam0 <- inits$lam0
    lamd <- lam0*exp(-D*D/(2*sigma*sigma))
    ll.y <- dpois(y.true2D,K*lamd*z,log=TRUE)
  }else if(obstype=="negbin"){
    lam0 <- inits$lam0
    theta <- inits$theta
    lamd <- lam0*exp(-D*D/(2*sigma*sigma))
    ll.y <- y.true2D*0
    for(i in 1:M){
      if(z[i]==1){
        ll.y[i,] <- dnbinom(y.true2D[i,],mu=lamd[i,],size=theta*K,log=TRUE)
      }
    }
  }else if(obstype=="hurdleZTPois"){
    p0 <- inits$p0
    lambda <- inits$lambda
    pd <- p0*exp(-D*D/(2*sigma*sigma))
    ll.y <- y.true*0
    for(i in 1:M){
      if(z[i]==1){
        for(j in 1:J){
          for(k in 1:K){
            if(y.true[i,j,k]==0){
              ll.y[i,j,k] <- log(1-pd[i,j])
            }else{
              ll.y[i,j,k] <- log(pd[i,j]) + log(dpois(y.true[i,j,k],lambda=lambda)/(1-exp(-lambda)))
            }
          }
        }
      }
    }
  }else{
    stop("obstype not recognized")
  }
  if(!is.finite(sum(ll.y)))stop("Starting obs model likelihood is not finite")
  if(obstype!="ramsey"){
    return(list(y.true2D=y.true2D,y.true3D=y.true,s=s,z=z,this.j=this.j,this.k=this.k,
                ID=ID,n.samples=n.samples,xlim=xlim,ylim=ylim))
  }else{
    return(list(y.true2D=y.true2D,y.true3D=y.true,s=s,z=z,
                jk.idx=jk.idx,xlim=xlim,ylim=ylim,ID=NA))
  }
}