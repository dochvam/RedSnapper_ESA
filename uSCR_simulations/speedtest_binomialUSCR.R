source("uSCR_simulations/SCR_sim_fn.R")

#### Default samplers ####

M_vec <- 1:4 * 20

for (i in 1:length(M_vec)) {
  for (t in 1:2) {
    test <- run_one_uSCR_simulation_binomial(iter = 1110 + i + t*10000, 
                                             M = M_vec[i], 
                                             niter = 50, nburnin = 0, 
                                             nchains = 1, thin = 1,
                                             prefix = "SPEEDTEST_upfiltered_")
  }
}



result_files <- list.files("intermediate/sim", pattern = "SPEEDTEST_simplebinomial", full.names = T)

result <- lapply(result_files, function(x) readRDS(x)$diagnostics) %>% 
  bind_rows() %>% 
  mutate(time_per_iter = as.numeric(mcmc_time) / 50) %>% 
  # bind_rows(validation_test$diagnostics %>% mutate(time_per_iter = as.numeric(mcmc_time) / 5)) %>% 
  mutate(M_sq = M^2) %>% 
  mutate(type = "Default")

result_uf_files <- list.files("intermediate/sim", pattern = "SPEEDTEST_upfiltered_simplebinomial", full.names = T)
result_uf <- lapply(result_uf_files, function(x) readRDS(x)$diagnostics) %>% 
  bind_rows() %>% 
  mutate(time_per_iter = as.numeric(mcmc_time) / 50) %>% 
  # bind_rows(validation_test$diagnostics %>% mutate(time_per_iter = as.numeric(mcmc_time) / 5)) %>% 
  mutate(M_sq = M^2) %>% 
  mutate(type = "Upfiltered")

ggplot(bind_rows(result, result_uf), aes(M_sq, as.numeric(time_per_iter), group = type, col = type)) +
  geom_point() +
  geom_smooth()

ggplot(result, aes(M, as.numeric(build_time))) +
  geom_point() +
  geom_smooth()

fit <- lm(as.numeric(time_per_iter) ~ M + M_sq, data = result)
summary(fit)

predict(fit, newdata = data.frame(M = c(100, 200, 600, 1000)) %>% mutate(M_sq = M^2))



# # Check: is our prediction (252 s per iter) matched if we bump M to 600?
# test <- run_one_uSCR_simulation_binomial(iter = 1110, 
#                                          M = 600, 
#                                          niter = 5, nburnin = 0, 
#                                          nchains = 1, thin = 1,
#                                          prefix = "SPEEDTEST_Check_simplebinomial")
# 
# validation_test <- readRDS("intermediate/sim/SPEEDTEST_Check_simplebinomial_uSCR_1110.RDS")
# validation_test$diagnostics
# 5 iterations should be ~1025 seconds


#### RW block samplers for each s[i,,] ####

# Let's try improving on this. What if we jointly sample s[i , , ] so we have to
# calculate the model fewer times?

# M_vec <- 1:8 * 20

for (i in 1:length(M_vec)) {
  for (t in 2:2) {
    test <- run_one_uSCR_simulation_binomial(iter = 1110 + i + t*10000, 
                                             M = M_vec[i], 
                                             niter = 50, nburnin = 0, 
                                             nchains = 1, thin = 1, 
                                             sampler_spec = "RW_block_1", 
                                             prefix = "SPEEDTEST_RWB_")
  }
}



result_files <- list.files("intermediate/sim", pattern = "SPEEDTEST_RWB_simplebinomial", full.names = T)

result2 <- lapply(result_files, function(x) readRDS(x)$diagnostics) %>% 
  bind_rows() %>% 
  mutate(time_per_iter = as.numeric(mcmc_time) / 50) %>% 
  mutate(M_sq = M^2) %>% 
  mutate(type = "RWB")

ggplot(bind_rows(result2, result), aes(group = type, col = type, 
                                       M^2, as.numeric(time_per_iter))) +
  geom_point() +
  geom_smooth()

ggplot(result, aes(M, as.numeric(build_time))) +
  geom_point() +
  geom_smooth()



fit <- lm(as.numeric(time_per_iter) ~ M + M_sq, data = result2)
summary(fit)

predict(fit, newdata = data.frame(M = c(100, 200, 600, 1000)) %>% mutate(M_sq = M^2))





#### RW block samplers for groups of 5 s[i,,] ####

# Let's try improving on this. What if we jointly sample s[i , , ] so we have to
# calculate the model fewer times?

# M_vec <- 1:5 * 20

for (i in 1:length(M_vec)) {
  for (t in 1:2) {
    test <- run_one_uSCR_simulation_binomial(iter = 1110 + i + t*10000, 
                                             M = M_vec[i], 
                                             niter = 50, nburnin = 0, 
                                             nchains = 1, thin = 1, 
                                             sampler_spec = "RW_block_2", 
                                             prefix = "SPEEDTEST_RWB2_upfiltered_")
  }
}



result_files <- list.files("intermediate/sim", pattern = "SPEEDTEST_RWB2_simplebinomial", full.names = T)
result3 <- lapply(result_files, function(x) readRDS(x)$diagnostics) %>% 
  bind_rows() %>% 
  mutate(time_per_iter = as.numeric(mcmc_time) / 50) %>% 
  # bind_rows(validation_test$diagnostics %>% mutate(time_per_iter = as.numeric(mcmc_time) / 10)) %>% 
  mutate(M_sq = M^2) %>% 
  mutate(type = "RWB2")

result_uf_files <- list.files("intermediate/sim", pattern = "SPEEDTEST_RWB2_upfiltered_simplebinomial", full.names = T)
result3_uf <- lapply(result_uf_files, function(x) readRDS(x)$diagnostics) %>% 
  bind_rows() %>% 
  mutate(time_per_iter = as.numeric(mcmc_time) / 50) %>% 
  # bind_rows(validation_test$diagnostics %>% mutate(time_per_iter = as.numeric(mcmc_time) / 10)) %>% 
  mutate(M_sq = M^2) %>% 
  mutate(type = "RWB2_Upfiltered")

ggplot(bind_rows(result2, result3, result), aes(group = type, col = type, 
                                       M^2, as.numeric(time_per_iter))) +
  geom_point() +
  geom_smooth()

ggplot(bind_rows(result, result_uf, result3, result3_uf), aes(group = type, col = type, 
                                       M^2, as.numeric(time_per_iter))) +
  geom_point() +
  geom_smooth()




fit <- lm(as.numeric(time_per_iter) ~ M + M_sq, data = result3)
summary(fit)

predict(fit, newdata = data.frame(M = c(100, 200, 600, 1000)) %>% mutate(M_sq = M^2))

#### Another idea: run the samplers for everything but "s" multiple times to improve mixing quality
# # Check: is our prediction (102 s per iter) matched if we bump M to 600?
test <- run_one_uSCR_simulation_binomial(iter = 7896,
                                         M = 600,
                                         niter = 10, nburnin = 0,
                                         nchains = 1, thin = 1, sampler_spec = "RW_block_2",
                                         prefix = "SPEEDTEST_Check_RWB2_")

validation_test <- readRDS("intermediate/sim/SPEEDTEST_Check_RWB2_simplebinomial_uSCR_7896.RDS")
validation_test$diagnostics







#### RW block samplers for groups of 10 s[i,,] ####

# Let's try improving on this. What if we jointly sample s[i , , ] so we have to
# calculate the model fewer times?

# M_vec <- 1:5 * 20

for (i in 1:length(M_vec)) {
  for (t in 1:2) {
    test <- run_one_uSCR_simulation_binomial(iter = 1110 + i + t*10000, 
                                             M = M_vec[i], 
                                             niter = 50, nburnin = 0, 
                                             nchains = 1, thin = 1, 
                                             sampler_spec = "RW_block_3", 
                                             prefix = "SPEEDTEST_RWB3_")
  }
}



result_files <- list.files("intermediate/sim", pattern = "SPEEDTEST_RWB3_simplebinomial", full.names = T)

summary <- lapply(result_files, function(x) readRDS(x)$summary) %>% 
  bind_rows()

result4 <- lapply(result_files, function(x) readRDS(x)$diagnostics) %>% 
  bind_rows() %>% 
  mutate(time_per_iter = as.numeric(mcmc_time) / 50) %>% 
  mutate(M_sq = M^2) %>% 
  mutate(type = "RWB3")

ggplot(bind_rows(result2, result3, result4, result), aes(group = type, col = type, 
                                                M^2, as.numeric(time_per_iter))) +
  geom_point() +
  geom_smooth()
