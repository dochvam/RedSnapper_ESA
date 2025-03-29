library(tidyverse)
library(nimble)


result_file <- "uSCR_real/joint_masked_posterior_NB.RDS"
distr <- "NB"
prefix_out <- "joint_masked_uSCR"

temp <- readRDS(result_file)


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
    stop("Ben, implement binomial")
  } else {
    stop("Unrecognized distribution")
  }
}

ESA_df <- data.frame(
  current_dir = c("Towards", "Away", "Perpendicular"),
  distr = distr,
  prefix = prefix_out,
  ESA_q50 = apply(esa_posterior, 2, median),
  ESA_q025 = apply(esa_posterior, 2, quantile, prob = 0.025),
  ESA_q975 = apply(esa_posterior, 2, quantile, prob = 0.975)
)

write_csv(ESA_df, paste0("ESA_estimates/", prefix_out, "_", distr, "_ESA.csv"))

