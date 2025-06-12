calc_ESA <- function(result_list, distr, prefix_out, esa_sigma_type, 
                     log_sigma_estimate = NULL, 
                     write.out = F, nframe = 41, thin = 10) {
  
  stopifnot(esa_sigma_type %in% c("20-min", 
                                  "NB attraction dist",
                                  "NB attraction dist plus 20-min",
                                  "as modeled"))
  
  attraction_dist_mean <- 24.3870801
  attraction_dist_sd <- 4.188949
  
  temp <- result_list
  
  samples_mtx <- as.matrix(temp$samples)
  
  resoln <- 5
  max_m <- 500
  cell_area <- resoln^2
  grid <- expand.grid(x = seq(-max_m, max_m, by = resoln),
                      y = seq(-max_m, max_m, by = resoln)) %>% 
    as.data.frame() %>% 
    mutate(dist = sqrt(x^2 + y^2))
  
  if (esa_sigma_type == "as modeled") {
    this_sigma <- samples_mtx[, "sigma"]
  } else if (esa_sigma_type == "NB attraction dist plus 20-min") {
    this_logsigma <- rnorm(nrow(samples_mtx), 
                           log_sigma_estimate$mean[log_sigma_estimate$t_interval == 20 & log_sigma_estimate$type == "mean"],
                           log_sigma_estimate$sd[log_sigma_estimate$t_interval == 20 & log_sigma_estimate$type == "mean"])
    
    attraction_dist <- rnorm(nrow(samples_mtx), attraction_dist_mean, attraction_dist_sd)
    
    this_sigma <- sqrt(exp(this_logsigma)^2 + attraction_dist^2)
  } else if (esa_sigma_type == "NB attraction dist") {
    attraction_dist <- rnorm(nrow(samples_mtx), attraction_dist_mean, attraction_dist_sd)
    this_sigma <- attraction_dist
  } else if (esa_sigma_type == "20-min") {
    this_logsigma <- rnorm(nrow(samples_mtx), 
                           log_sigma_estimate$mean[log_sigma_estimate$t_interval == 20 & log_sigma_estimate$type == "mean"],
                           log_sigma_estimate$sd[log_sigma_estimate$t_interval == 20 & log_sigma_estimate$type == "mean"])
    this_sigma <- exp(this_logsigma)
  }
  
  
  # Matrix of posterior distributions for the ESA
  # Each column represents a current direction
  esa_posterior <- matrix(NA, nrow = nrow(samples_mtx), ncol = 3)
  
  # Loop over the posterior samples
  for (i in 10 * 1:(nrow(samples_mtx)/thin)) {
    # What's the probability that we get at least one detection at a trap placed 
    # at 0,0 given the fish's centroid is in each grid coord?
    if (distr == "Pois") {
      stop()
      for (j in 1:3) {
        lambda <- samples_mtx[i, paste0("lam0[", j, "]")] * exp(-grid$dist^2 / (2 * samples_mtx[i, "sigma"]^2))
        prob_one_det <- 1-exp(-lambda)
        esa_posterior[i, j] <- sum(prob_one_det * cell_area)
      }
      
    } else if (distr == "NB") {
      for (j in 1:3) {
        stop()
        lambda <- samples_mtx[i, paste0("lam0[", j, "]")] * exp(-grid$dist^2 / (2 * samples_mtx[i, "sigma"]^2))
        size_NB <- 1 / samples_mtx[i, "phi_CT"]
        prob_NB <- 1 / (1 + samples_mtx[i, "phi_CT"] * lambda)
        prob_one_det <- 1 - dnbinom(0, size = size_NB, prob_NB)
        esa_posterior[i, j] <- sum(prob_one_det * cell_area)
      }
    } else if (distr == "Binom") {
      for (j in 1:3) {

        p <- samples_mtx[i, paste0("p0[", j, "]")] * exp(-grid$dist^2 /(2*this_sigma[i]^2))
        p_k <- (1 - (1-p)^nframe)
        esa_posterior[i, j] <- sum(p_k * cell_area)
        
        # esa_posterior[i, j] <- (samples_mtx[i, "sigma"] * sqrt(2 * samples_mtx[i, paste0("p0[", j, "]")]))^2 * base::pi # from half-normal EDD
      }
    } else {
      stop("Unrecognized distribution")
    }
  }
  
  ESA_df <- data.frame(
    current_dir = c("Towards", "Away", "Across"),
    distr = distr,
    mean = apply(esa_posterior, 2, mean, na.rm = T),
    sd = apply(esa_posterior, 2, sd, na.rm = T),
    ESA_q50 = apply(esa_posterior, 2, median, na.rm = T),
    ESA_q025 = apply(esa_posterior, 2, quantile, prob = 0.025, na.rm = T),
    ESA_q975 = apply(esa_posterior, 2, quantile, prob = 0.975, na.rm = T)
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
