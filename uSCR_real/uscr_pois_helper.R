
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
                         spatial_beta,
                         habitatMask) {
  s <- array(NA, dim =c(n, 3, 2))
  for (i in 1:n) {
    for (ti in 1:t) {
      s[i, ti, ] <- rHabDistr_asCovar(1, xmax = xmax, xmin = xmin,
                                      ymax = ymax, ymin = ymin, resoln = resoln,
                                      phi = exp(spatial_beta * log(habitatMask)))
      
    }
  }
  s
}