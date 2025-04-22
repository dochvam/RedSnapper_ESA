
#### nimbleFunction GetDetectionRate ####
GetDetectionRate <- nimbleFunction(
  run = function(s = double(1), lam0=double(1), 
                 sigma=double(0), 
                 current_dir = double(1),
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      ans <- lam0[current_dir]*exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)


dHabDistr_asCovar <- nimbleFunction(run = function(x = double(1),
                                                   xmax = double(0),
                                                   xmin = double(0),
                                                   ymax = double(0),
                                                   ymin = double(0),
                                                   resoln = double(0),
                                                   phi = double(2),
                                                   # beta = double(0),
                                                   log = logical(0)) {
  if (x[1] < xmin) return(-Inf)
  if (x[1] > xmax) return(-Inf)
  if (x[2] < ymin) return(-Inf)
  if (x[2] > ymax) return(-Inf)
  
  # Convert s to row/col
  row <- trunc((x[1] - xmin) / resoln) + 1
  col <- trunc((x[2] - ymin) / resoln) + 1
  
  # phi <- exp(beta * habitatMask) # hm is log-scale covariate
  
  if (log) return(log(phi[row, col] / sum(phi)))
  else return(phi[row, col] / sum(phi))
  returnType(double(0))
}, buildDerivs = FALSE)
# ctest <- compileNimble(dHabDistr_asCovar)



rHabDistr_asCovar <- nimbleFunction(run = function(n = double(0),
                                                   xmax = double(0),
                                                   xmin = double(0),
                                                   ymax = double(0),
                                                   ymin = double(0),
                                                   resoln = double(0),
                                                   phi = double(2)#,
                                                   # beta = double(0)
){
  returnType(double(1))
  
  # phi <- exp(beta * log(habitatMask)) # log(hm) is covariate
  prob_mtx <- phi / sum(phi)
  
  # Randomly select a row based on their relative weights
  rowProbs <- numeric(dim(prob_mtx)[1])
  for (i in 1:(dim(prob_mtx)[1])) {
    rowProbs[i] <- sum(prob_mtx[i, ])
  }
  row <- rcat(1, prob = rowProbs)
  
  # Randomly select a column from that row
  col <- rcat(1, prob = prob_mtx[row, ])
  
  # # Generate a random coord. uniformly within that cell
  x <- c(
    xmin + resoln * (row - 1) + runif(1, 0, resoln),
    ymin + resoln * (col - 1) + runif(1, 0, resoln)
  )
  return(x)
}, buildDerivs = FALSE)


calcPctFishInROVbuffer <- nimbleFunction(run = function(phi = double(2), 
                                                        weights = double(1), 
                                                        rov_cell_xvec = double(1),
                                                        rov_cell_yvec = double(1),
                                                        n = integer(0)) {
  sum <- 0
  for (i in 1:n) {
    sum <- sum + weights[i] * phi[rov_cell_xvec[i], rov_cell_yvec[i]]
  }
  return(sum)
  returnType(double())
})



initialize_s <- function(n, t,
                         xmax,
                         xmin,
                         ymax,
                         ymin,
                         resoln,
                         habitatMask,
                         spatial_beta) {
  s <- array(NA, dim =c(n, 3, 2))
  for (i in 1:n) {
    for (ti in 1:t) {
      # browser()
      s[i, ti, ] <- rHabDistr_asCovar(1, xmax = xmax, xmin = xmin,
                                      ymax = ymax, ymin = ymin, resoln = resoln,
                                      phi = exp(spatial_beta * log(habitatMask)))
      
    }
  }
  s
}



dPoisBinom <- nimbleFunction(run = function(x = double(0),
                                            p = double(1),
                                            log = logical(0)) {
  # p <- p[p > 1e-10] # Ignore trials for which p < 1e-10
  # if (length(p) == 0) {
  #   if (x == 0) {
  #     prob <- 1
  #   } else {
  #     prob <- 0
  #   }
  #   if (log) return(log(prob))
  #   else return(prob)
  # }
  
  n <- length(p)
  if (x > n) {
    prob <- 0
  } else {
    prob_vec <- numeric(n + 1)
    prob_vec[1] <- 1
    
    for (i in 1:n) {
      q <- 1 - p[i]
      prob_vec[2:(i + 1)] <- prob_vec[2:(i + 1)] * q + prob_vec[1:i] * p[i]
      prob_vec[1] <- prob_vec[1] * q
    }
    prob <- prob_vec[x + 1]
  }
  
  
  if (log) return(log(prob))
  else return(prob)
  returnType(double(0))
})

rPoisBinom <- nimbleFunction(run = function(n = double(0),
                                            p = double(1)) {
  if (n != 1) stop("n must equal 1 in rPoisBinom")
  return(sum(rbinom(n = length(p), size = 1, prob = p)))
  returnType(double(0))
})

dPoisBinom_wReps <- nimbleFunction(run = function(x = double(0),
                                                  p = double(1),
                                                  reps = integer(0),
                                                  log = logical(0)) {
  # 
  # p_filtered <- p
  p_filtered <- p[p > 1e-5] # Ignore trials for which p < 1e-10
  if (length(p_filtered) == 0) {
    if (x == 0) {
      prob <- 1
    } else {
      prob <- 0
    }
    if (log) return(log(prob))
    else return(prob)
  }
  
  max_successes <- reps * length(p_filtered)
  probs <- numeric(max_successes + 1)
  probs[1] <- 1
  
  if (x > max_successes) {
    prob <- 0
  } else {
    
    for (i in 1:length(p_filtered)) {
      q <- 1 - p_filtered[i]
      
      new_probs <- numeric(max_successes + 1)
      
      for (k in 0:reps) {
        binom_prob <- dbinom(k, reps, p_filtered[i])
        new_probs[(k+1):(max_successes+1)] <- 
          new_probs[(k+1):(max_successes+1)] + probs[1:(max_successes + 1 - k)] * binom_prob
      }
      
      probs <- new_probs
    }
    prob <- probs[x + 1]
  }
  
  
  if (log) return(log(prob))
  else return(prob)
  returnType(double(0))
})


rPoisBinom_wReps <- nimbleFunction(run = function(n = double(0),
                                                  p = double(1),
                                                  reps = integer(0)) {
  if (n != 1) stop("n must equal 1 in rPoisBinom")
  return(sum(rbinom(n = length(p), size = reps, prob = p)))
  returnType(double(0))
})



##### Samplers #####
# Custom RJ sampler from Sally Paganin + Daniel Eacker + Perry de Valpine
my_sampler_MV_RJ_indicator <- nimbleFunction(
  name = 'my_sampler_MV_RJ_indicator',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## note: target is the indicator variable,
    ## control$targetNode is the variable conditionally in the model
    ## control list extraction
    coefNode      <- model$expandNodeNames(control$targetNode, returnScalarComponents = TRUE)
    # coefNodes     <- control$targetNode
    
    proposalScale <- control$scale
    proposalMean  <- control$mean
    len_coefNode <- length(coefNode) # It is better to do this in setup code and use below
    
    
    ## node list generation
    calcNodes <- model$getDependencies(c(coefNode, target))
    calcNodesReduced <- model$getDependencies(target)
  },
  run = function() {
    currentIndicator <- model[[target]]
    if(currentIndicator == 0) {   ## propose addition of coefNode
      currentLogProb <- model$getLogProb(calcNodesReduced)
      proposalCoef <- numeric(len_coefNode)
      logProbForwardProposal <- 0
      for(l in 1:len_coefNode) {
        
        proposalCoef[l] <- rnorm(1, proposalMean, proposalScale)
        logProbForwardProposal <- logProbForwardProposal + dnorm(proposalCoef[l], proposalMean, proposalScale, log = TRUE)
      }
      values(model, coefNode) <<- proposalCoef
      
      model[[target]] <<- 1
      proposalLogProb <- model$calculate(calcNodes)
      logAcceptanceProb <- proposalLogProb - currentLogProb - logProbForwardProposal
    } else {                      ## propose removal of coefNode
      currentLogProb <- model$getLogProb(calcNodes)
      currentCoef <-  values(model, coefNode)
      logProbReverseProposal<- 0
      for(l in 1:len_coefNode) {
        
        logProbReverseProposal <- logProbReverseProposal + dnorm(currentCoef[l], proposalMean, proposalScale, log = TRUE)
      }
      values(model, coefNode) <<- rep(0, len_coefNode)
      
      model[[target]] <<- 0
      model$calculate(calcNodes)
      logAcceptanceProb <- model$getLogProb(calcNodesReduced) - currentLogProb + logProbReverseProposal
    }
    accept <- decide(logAcceptanceProb)
    if(accept) { copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    } else     { copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE) }
  },
  methods = list(
    reset = function() { }
  )
)

my_sampler_RJ_toggled <- nimbleFunction(
  name = 'sampler_RJ_toggled',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## control list extraction
    fixedValue  <- if(!is.null(control$fixedValue))  control$fixedValue  else 0
    samplerType <- if(!is.null(control$samplerType)) control$samplerType else stop('must provide \'samplerType\' control list element to RJ_toggled sampler')
    ## nested function and function list definitions
    samplerFL <- nimbleFunctionList(sampler_BASE)
    samplerFL[[1]] <- samplerType$buildSampler(model = model, mvSaved = mvSaved)
  },
  run = function() {
    if(!all(model[[target]] != fixedValue))
      samplerFL[[1]]$run()
  },
  methods = list(
    reset = function() {
      samplerFL[[1]]$reset()
    }
  )
)


