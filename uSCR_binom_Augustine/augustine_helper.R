library(nimble)


#------------------------------------------------------------------
# Everything in this script was written by Ben Augustine
#------------------------------------------------------------------
dBernoulliVector <- nimbleFunction(
  run = function(x = double(1), pd = double(1), K1D=double(1), z = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      logProb <- sum(dbinom(x, prob = pd, size = K1D, log = TRUE))
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rBernoulliVector <- nimbleFunction(
  run = function(n = integer(0),pd = double(1),K1D=double(1), z = double(0)) {
    returnType(double(1))
    J <- nimDim(pd)[1]
    out <- numeric(J,value=0)
    return(out)
  }
)

Getcapcounts <- nimbleFunction(
  run = function(ID=double(1),M=double(0)){
    returnType(double(1))
    n.samples <- nimDim(ID)[1]
    capcounts <- numeric(M, value = 0)
    for(l in 1:n.samples){
      capcounts[ID[l]] <- capcounts[ID[l]] + 1
    }
    return(capcounts)
  }
)

Getncap <- nimbleFunction(
  run = function(capcounts=double(1)){
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

## sampler to update y[1:M,1:J,1:K] (3D data), then fed back to nimble as y[1:M,1:J] (2D data)
IDSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    this.j <- control$this.j
    this.k <- control$this.k
    n.samples <- control$n.samples
    M <- control$M
    J <- control$J
    K <- control$K
    IDups <- control$IDups
    calcNodes <- model$getDependencies(c("y.true","ID"))
  },
  
  run = function() {
    z <- model$z
    # y.true <- model$y.true
    ID.curr <- model$ID
    
    #build y.true in 3D
    y.true <- array(0,dim=c(M,J,K))
    for(l in 1:n.samples){
      y.true[ID.curr[l],this.j[l],this.k[l]] <- 1
    }
    
    #precalculate log likelihoods 3D
    ll.y <- array(0,dim=c(M,J,K))
    for(i in 1:M){
      if(z[i]==1){
        for(j in 1:J){
          for(k in 1:K){
            ll.y[i,j,k] <-  dbinom(y.true[i,j,k],size=1,prob=model$pd[i,j], log = TRUE)
          }
        }
      }
    }
    ll.y.cand <- ll.y
    ID.cand <- ID.curr
    y.true.cand <- y.true
    
    for(up in 1:IDups){
      for(l in 1:n.samples){ #loop over samples
        #proposal distribution
        propprobs <- model$pd[1:M,this.j[l]]*z[1:M]
        sumpropprobs <- sum(propprobs)
        propprobs <- propprobs/sumpropprobs
        ID.cand[l] <- rcat(1,prob=propprobs)
        if(ID.cand[l]!=ID.curr[l]){
          swapped <- c(ID.curr[l],ID.cand[l])
          #new sample proposal probabilities
          forprob <- propprobs[swapped[2]]
          backprob <- propprobs[swapped[1]]
          
          #does cand guy have a sample at this j-k already?
          occupied <- sum(ID.curr==ID.cand[l]&this.j==this.j[l]&this.k==this.k[l])==1
          if(occupied){ #swap samples, y.true and likelihood don't change
            #need to exchange this sample with focal guy
            swap2 <- which(ID.curr==ID.cand[l]&this.j==this.j[l]&this.k==this.k[l])
            ID.cand[swap2] <- ID.curr[l]
          }else{
            #new y.true's - move sample from ID to ID.cand
            y.true.cand[ID.curr[l],this.j[l],this.k[l]] <- y.true[ID.curr[l],this.j[l],this.k[l]] - 1
            y.true.cand[ID.cand[l],this.j[l],this.k[l]] <- y.true[ID.cand[l],this.j[l],this.k[l]] + 1
          }
          
          ll.y.cand[swapped[1],this.j[l],this.k[l]] <- dbinom(y.true.cand[swapped[1],this.j[l],this.k[l]],size=1,prob=model$pd[swapped[1],this.j[l]],log=TRUE)
          ll.y.cand[swapped[2],this.j[l],this.k[l]] <- dbinom(y.true.cand[swapped[2],this.j[l],this.k[l]],size=1,prob=model$pd[swapped[2],this.j[l]],log=TRUE)
          
          #select sample to move proposal probabilities
          focalprob <- sum(ID.curr==swapped[1]&this.j==this.j[l]&this.k==this.k[l])/n.samples
          focalbackprob <- sum(ID.cand==swapped[2]&this.j==this.j[l]&this.k==this.k[l])/n.samples
          
          #sum log likelihoods and do MH step
          lp_initial <- ll.y[swapped[1],this.j[l],this.k[l]] + ll.y[swapped[2],this.j[l],this.k[l]]
          lp_proposed <- ll.y.cand[swapped[1],this.j[l],this.k[l]] + ll.y.cand[swapped[2],this.j[l],this.k[l]]
          log_MH_ratio <- (lp_proposed+log(backprob)+log(focalbackprob)) - (lp_initial+log(forprob)+log(focalprob))
          
          accept <- decide(log_MH_ratio)
          if(accept){
            y.true[swapped[1],this.j[l],this.k[l]] <- y.true.cand[swapped[1],this.j[l],this.k[l]]
            y.true[swapped[2],this.j[l],this.k[l]] <- y.true.cand[swapped[2],this.j[l],this.k[l]]
            ll.y[swapped[1],this.j[l],this.k[l]] <- ll.y.cand[swapped[1],this.j[l],this.k[l]]
            ll.y[swapped[2],this.j[l],this.k[l]] <- ll.y.cand[swapped[2],this.j[l],this.k[l]]
            ID.curr[l] <- ID.cand[l]
            if(occupied){
              ID.curr[swap2] <- ID.cand[swap2]
            }
          }else{ #set back
            y.true.cand[swapped[1],this.j[l],this.k[l]] <- y.true[swapped[1],this.j[l],this.k[l]]
            y.true.cand[swapped[2],this.j[l],this.k[l]] <- y.true[swapped[2],this.j[l],this.k[l]]
            ll.y.cand[swapped[1],this.j[l],this.k[l]] <- ll.y[swapped[1],this.j[l],this.k[l]]
            ll.y.cand[swapped[2],this.j[l],this.k[l]] <- ll.y[swapped[2],this.j[l],this.k[l]]
            ID.cand[l] <- ID.curr[l]
            if(occupied){
              ID.cand[swap2] <- ID.curr[swap2]
            }
          }
        }else{ #set back
          ID.cand[l] <- ID.curr[l]
        }
      }
    }
    
    #rebuild y.true in 2D
    y.true2D <- matrix(0,M,J)
    for(l in 1:n.samples){
      y.true2D[ID.curr[l],this.j[l]] <- y.true2D[ID.curr[l],this.j[l]] + 1
    }
    
    #put everything back into the model$stuff
    model$y.true <<- y.true2D
    model$ID <<- ID.curr
    model$calculate(calcNodes) #update logprob
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)

#Required custom update for N/z
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    M <- control$M
    z.ups <- control$z.ups
    y.nodes <- control$y.nodes
    pd.nodes <- control$pd.nodes
    N.node <- control$N.node
    z.nodes <- control$z.nodes
    calcNodes <- control$calcNodes
  },
  run = function() {
    for(up in 1:z.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected call
      if(updown==0){#subtract
        # find all z's currently on
        z.on <- which(model$z==1)
        n.z.on <- length(z.on)
        pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
        pick <- z.on[pick]
        #prereject turning off individuals currently allocated samples
        if(model$capcounts[pick]>0){#is this an individual with samples?
          reject <- TRUE
        }
        if(!reject){
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick])
          
          #propose new N/z
          model$N[1] <<-  model$N[1] - 1
          model$z[pick] <<- 0
          
          #turn pd off
          model$calculate(pd.nodes[pick])
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick]) #will always be 0
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          accept <- decide(log_MH_ratio)
          
          if(accept) {
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["pd",1][pick,] <<- model[["pd"]][pick,]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["pd"]][pick,] <<- mvSaved["pd",1][pick,]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }else{#add
        if(model$N[1] < M){ #cannot update if z maxed out. Need to raise M
          z.off <- which(model$z==0)
          n.z.off <- length(z.off)
          pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
          pick <- z.off[pick]
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick]) #will always be 0
          
          #propose new N/z
          model$N[1] <<-  model$N[1] + 1
          model$z[pick] <<- 1
          
          #turn pd on
          model$calculate(pd.nodes[pick])
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick])
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["pd",1][pick,] <<- model[["pd"]][pick,]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["pd"]][pick,] <<- mvSaved["pd",1][pick,]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }
    }
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)


sSampler <- nimbleFunction(
  # name = 'sampler_RW',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    i<-control$i    
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
      model$s[i, 1:2] <<- c(runif(1, xlim[1], xlim[2]), runif(1, ylim[1], ylim[2]))
      model$calculate(calcNodes)
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    }else{#MH
      s.cand=c(rnorm(1,model$s[i,1],scale), rnorm(1,model$s[i,2],scale))
      inbox= s.cand[1]< xlim[2] & s.cand[1]> xlim[1] & s.cand[2] < ylim[2] & s.cand[2] > ylim[1]
      if(inbox){
        model_lp_initial <- model$getLogProb(calcNodes)
        model$s[i, 1:2] <<- s.cand
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