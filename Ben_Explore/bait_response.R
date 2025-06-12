library(tidyverse)
library(mgcv)
library(nimble)
library(MCMCvis)

responses <- read_csv("Data/GAM_input_rs_30min.csv")

bait_response_model <- nimbleCode({
  
  for (i in 1:nobs) {
    response_sigma_indiv[i] <- exp(log(response_sigma) + individual_ranef[fish[i]])
    p[i] <- exp(-initial_dist[i]^2 / (2 * response_sigma_indiv[i]^2))
    y[i] ~ dbern(p[i])
  }
  
  response_sigma ~ dunif(0, 500)
  ranef_sd ~ dgamma(1, 1)
  for (i in 1:nfish) {
    individual_ranef[i] ~ dnorm(0, sd = ranef_sd)
  }
})


mod <- nimbleModel(bait_response_model,
                   constants = list(nfish = length(unique(responses$num)),
                                    nobs = nrow(responses),
                                    initial_dist = responses$initial.dist,
                                    fish = as.numeric(as.factor(responses$num))),
                   data = list(y = responses$response30),
                   inits = list(response_sigma = 50, ranef_sd = 0.1,
                                individual_ranef = rnorm(length(unique(responses$num)), 0, 0.1)))

samples <- nimbleMCMC(mod, samplesAsCodaMCMC = TRUE, 
                      niter = 10000, nburnin = 5000, nchains = 2)
summary <- MCMCsummary(samples)



predictions <- read_csv("Data/GAM_RS_results_30m_response.csv") %>% 
  mutate(type = "GAM")

predictions_hn <- predictions[, "initial.dist"] %>% 
  mutate(type = "Half-normal")
predictions_hn$response.30m <- NA
predictions_hn$upperCI <- NA
predictions_hn$lowerCI <- NA

sigma_samples <- as.numeric(as.matrix(samples[, "response_sigma"]))

# Predict at each distance
for (i in 1:nrow(predictions_hn)) {
  p_vec <- exp(-predictions_hn$initial.dist[i]^2 / (2 * sigma_samples^2))
  predictions_hn$response.30m[i] <- quantile(p_vec, probs = 0.5)
  predictions_hn$upperCI[i] <- quantile(p_vec, probs = 0.975)
  predictions_hn$lowerCI[i] <- quantile(p_vec, probs = 0.025)
}


bind_rows(predictions, predictions_hn) %>% 
  ggplot(aes(initial.dist, response.30m, ymin = lowerCI, ymax = upperCI,
             group = type, fill = type, col = type)) +
  geom_ribbon(alpha = 0.2, col = NA) + 
  geom_line(linewidth = 0.7) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  xlim(0, 150) +
  scale_fill_manual(values = c("black", "darkred")) +
  scale_color_manual(values = c("black", "darkred"))



