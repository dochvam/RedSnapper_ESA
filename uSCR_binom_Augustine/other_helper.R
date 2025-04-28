library(nimble)

GetDetectionProb_wDates <- nimbleFunction(
  run = function(s = double(1), p0=double(1), sigma=double(0), 
                 datevec = double(1), idate = double(0),
                 X_mtx=double(2), current_dir = double(1),
                 J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
      d2 <- numeric(J)
      ans <- numeric(J)
      for (j in 1:J) {
        if (datevec[j] == idate) {
          d2[j] <- (s[1]-X_mtx[j,1])^2 + (s[2]-X_mtx[j,2])^2
          ans[j] <- p0[current_dir[j]] * exp(-d2[j]/(2*sigma^2))
        } else {
          d2[j] <- 0
          ans[j] <- 0
        }
      }
      return(ans)
    }
  }
)


initialize_s <- function(M,
                         xlim,
                         ylim,
                         resoln,
                         habitatMask,
                         spatial_beta) {
  s <- array(NA, dim =c(M, 2))
  for (i in 1:M) {
    s[i, ] <- rHabDistr_asCovar(1, xmax = xlim[2], xmin = xlim[1],
                                ymax = ylim[2], ymin = ylim[1], resoln = resoln,
                                phi = exp(spatial_beta * log(habitatMask)))
  }
  s
}



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



# sim_counts <- function(M, psi, p0, sigma, X, hm, spatial_beta) {
sim_counts <- function(M, psi, p0, sigma, current_dir, X_mtx, X_datevec, xlim, ylim, J, K) {
  
  N <- rbinom(1, M / 3, psi)
  
  idate <- rep(1:3, each = N)
  
  s <- matrix(NA, nrow = N*3, ncol = 2)
  s[,1] <- runif(N*3, xlim[1], xlim[2])
  s[,2] <- runif(N*3, ylim[1], ylim[2])
  
  y <- matrix(0, nrow(X_mtx), K)
  
  for (i in 1:nrow(s)) {
    # this_inds <- which(X_datevec == i_datevec[i])
    this_p <- GetDetectionProb_wDates(s = s[i,1:2], X_mtx = X_mtx, 
                                      datevec = X_datevec, idate = idate[i],
                                      current_dir = current_dir,
                                      J=nrow(X_mtx), sigma=sigma, p0=p0, z=1)
    for (k in 1:K) {
      y[, k] <- y[, k] + rbinom(J, 1, this_p)
    }
  }
  
  n.samples <- sum(y)
  
  this.j <- this.k <- rep(NA,n.samples)
  idx <- 1
  for(j in 1:J){ #then traps
    for(k in 1:K){ #then occasions
      if(y[j,k]>0){ #is there at least one sample here?
        for(l in 1:y[j,k]){ #then samples
          this.j[idx] <- j
          this.k[idx] <- k
          idx <- idx+1
        }
      }
    }
  }
  
  
  return(list(
    y = y, this.j = this.j, this.k = this.k, n.samples = n.samples,
    trueN = N, trueN3 = N*3
  ))
}



sim_counts_wHabCovar <- function(M, psi, p0, sigma, X_mtx, X_datevec, 
                                 current_dir,
                                 xlim, ylim, J, K, resoln, habitatMask,
                                 spatial_beta) {
  
  N <- rbinom(1, M / 3, psi)
  
  idate <- rep(1:3, each = N)
  
  s <- initialize_s(
    M = N*3,
    xlim = xlim,
    ylim = ylim,
    resoln = resoln,
    habitatMask = habitatMask,
    spatial_beta = spatial_beta
  )
  
  y <- matrix(0, nrow(X_mtx), K)
  
  for (i in 1:nrow(s)) {
    # this_inds <- which(X_datevec == i_datevec[i])
    this_p <- GetDetectionProb_wDates(s = s[i,1:2], X_mtx = X_mtx, 
                                      datevec = X_datevec, idate = idate[i],
                                      current_dir = current_dir,
                                      J=nrow(X_mtx), sigma=sigma, p0=p0, z=1)
    for (k in 1:K) {
      y[, k] <- y[, k] + rbinom(J, 1, this_p)
    }
  }
  
  n.samples <- sum(y)
  
  this.j <- this.k <- rep(NA,n.samples)
  idx <- 1
  for(j in 1:J){ #then traps
    for(k in 1:K){ #then occasions
      if(y[j,k]>0){ #is there at least one sample here?
        for(l in 1:y[j,k]){ #then samples
          this.j[idx] <- j
          this.k[idx] <- k
          idx <- idx+1
        }
      }
    }
  }
  
  
  return(list(
    y = y, this.j = this.j, this.k = this.k, n.samples = n.samples,
    trueN = N, trueN3 = N*3, true_s = s
  ))
}

initialize_ytrue <- function(M, z_init, s_init, this.j, this.k, X_mtx, current_dir,
                             X_datevec, idate, n.samples, sigma_init, p0_init, K) {
  # ytrue2D needs to be a matrix of [M] x [J], compatible with z, with sum = n.samples
  
  detprobs <- matrix(0, length(idate), nrow(X_mtx))
  for (i in 1:nrow(detprobs)) {
    detprobs[i,] <- GetDetectionProb_wDates(s = s_init[i,1:2], X_mtx = X_mtx, 
                                           datevec = X_datevec, idate = idate[i],
                                           current_dir = current_dir,
                                           J=nrow(X_mtx), sigma=sigma_init, p0=p0_init, z=z_init[i])
  }
  
  y_true_3D <- array(0, dim = c(M, nrow(X_mtx), K))
  ID <- numeric(length(this.j))
  
  for (o in 1:length(this.j)) {
    probs <- detprobs[, this.j[o]] / sum(detprobs[, this.j[o]])
    impossible <- (y_true_3D[, this.j[o], this.k[o]] > 0)
    probs[impossible] <- 0
    probs <- probs / sum(probs)
    
    ind <- sample(1:(M), size = 1, prob = probs)
    y_true_3D[ind, this.j[o], this.k[o]] <- y_true_3D[ind, this.j[o], this.k[o]] + 1
    
    ID[o] <- ind
  }
  
  return(list(
    ytrue3D = y_true_3D,
    ytrue2D = apply(y_true_3D, c(1,2), sum),
    ID = ID
  ))
}


sim_ROV <- function(nROV,
                    rb_weights,
                    rov_cell_xvec,
                    rov_cell_yvec,
                    rbe, 
                    rbs,
                    spatial_beta,
                    habitatMask,
                    sim.out) {
  
  ROV_obs <- numeric(nROV)
  phi <- exp(spatial_beta * log(habitatMask))
  
  for (i in 1:nROV) {
    pctFishInROVbuffer <- calcPctFishInROVbuffer(phi = phi, 
                                                 weights = rb_weights[rbs[i]:rbe[i]], 
                                                 rov_cell_xvec = rov_cell_xvec[rbs[i]:rbe[i]],
                                                 rov_cell_yvec = rov_cell_yvec[rbs[i]:rbe[i]],
                                                 n = rbe[i] - rbs[i] + 1)
    
    
    ROV_obs[i] <- rpois(1, pctFishInROVbuffer * sim.out$trueN)
  }
  
  ROV_obs
}




