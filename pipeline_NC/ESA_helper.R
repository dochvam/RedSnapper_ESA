calc_ESA <- function(result_list, distr, prefix_out, write.out = F) {
  
  temp <- result_list
  
  samples_mtx <- as.matrix(temp$samples)
  
  resoln <- 5
  max_m <- 500
  cell_area <- resoln^2
  grid <- expand.grid(x = seq(-max_m, max_m, by = resoln),
                      y = seq(-max_m, max_m, by = resoln)) %>% 
    as.data.frame() %>% 
    mutate(dist = sqrt(x^2 + y^2))
  
  # Matrix of posterior distributions for the ESA
  # Each column represents a current direction
  esa_posterior <- matrix(NA, nrow = nrow(samples_mtx), ncol = 3)
  
  # Loop over the posterior samples
  for (i in 1:nrow(samples_mtx)) {
    # What's the probability that we get at least one detection at a trap placed 
    # at 0,0 given the fish's centroid is in each grid coord?
    if (distr == "Pois") {
      for (j in 1:3) {
        lambda <- samples_mtx[i, paste0("lam0[", j, "]")] * exp(-grid$dist^2 / (2 * samples_mtx[i, "sigma"]^2))
        prob_one_det <- 1-exp(-lambda)
        esa_posterior[i, j] <- sum(prob_one_det * cell_area)
      }
      
    } else if (distr == "NB") {
      for (j in 1:3) {
        lambda <- samples_mtx[i, paste0("lam0[", j, "]")] * exp(-grid$dist^2 / (2 * samples_mtx[i, "sigma"]^2))
        size_NB <- 1 / samples_mtx[i, "phi_CT"]
        prob_NB <- 1 / (1 + samples_mtx[i, "phi_CT"] * lambda)
        prob_one_det <- 1 - dnbinom(0, size = size_NB, prob_NB)
        esa_posterior[i, j] <- sum(prob_one_det * cell_area)
      }
    } else if (distr == "Binom") {
      for (j in 1:3) {
        esa_posterior[i, j] <- (samples_mtx[i, "sigma"] * sqrt(2 * samples_mtx[i, paste0("p0[", j, "]")]))^2 * base::pi # from half-normal EDD
      }
    } else {
      stop("Unrecognized distribution")
    }
  }
  
  ESA_df <- data.frame(
    current_dir = c("Towards", "Away", "Across"),
    distr = distr,
    ESA_q50 = apply(esa_posterior, 2, median),
    ESA_q025 = apply(esa_posterior, 2, quantile, prob = 0.025),
    ESA_q975 = apply(esa_posterior, 2, quantile, prob = 0.975)
  )
  
  if (write.out) {
    write_csv(ESA_df, paste0("ESA_estimates/", prefix_out, "_", distr, "_ESA.csv"))
  }
  
  ESA_df
}


calc_95UA <- function(result_file, distr, prefix_out) {
  temp <- readRDS(result_file)
  
  sigma_samples <- as.numeric(unlist(temp$samples[, "sigma"]))
  
  ua_samples <- base::pi * (sigma_samples * 2.447)^2
  
  ua_df <- data.frame(type = "95% UA", distr = distr, prefix = prefix_out,
                      UA_q50  = median(ua_samples),
                      UA_q025 = quantile(ua_samples, prob = 0.025),
                      UA_q975 = quantile(ua_samples, prob = 0.975)
  )
  
  ua_df
}


get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
