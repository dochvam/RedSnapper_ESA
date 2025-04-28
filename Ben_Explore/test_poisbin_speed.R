source("uSCR_real/uscr_binom_helper.R")
cPoisBinom <- compileNimble(dPoisBinom)

# Test various dpoisbin implementations

pp_mtx <- matrix(runif(500 * 5000, 0, 1), nrow = 5000)


test_poisbinom <- numeric(nrow(pp_mtx))

test_PoissonBinomial_DFFT <- numeric(nrow(pp_mtx))
test_PoissonBinomial_cv <- numeric(nrow(pp_mtx))
test_PoissonBinomial_rec <- numeric(nrow(pp_mtx))
test_PoissonBinomial_char <- numeric(nrow(pp_mtx))

# test_NIMBLE_uc <- numeric(nrow(pp_mtx))
test_NIMBLE_c <- numeric(nrow(pp_mtx))



poisbinom_time <- system.time(
  for (i in 1:nrow(pp_mtx)) {
    test_poisbinom[i] <- poisbinom::dpoisbinom(250, pp_mtx[i, ], log_d = TRUE)
  }
)

PoissonBinomial_cv_time <- system.time(
  for (i in 1:nrow(pp_mtx)) {
    test_PoissonBinomial_cv[i] <- PoissonBinomial::dpbinom(250, pp_mtx[i, ], method = "Convolve", log = T)
  }
)

PoissonBinomial_DFFT_time <- system.time(
  for (i in 1:nrow(pp_mtx)) {
    test_PoissonBinomial_DFFT[i] <- PoissonBinomial::dpbinom(250, pp_mtx[i, ], method = "DivideFFT", log = T)
  }
)

PoissonBinomial_char_time <- system.time(
  for (i in 1:nrow(pp_mtx)) {
    test_PoissonBinomial_char[i] <- PoissonBinomial::dpbinom(x = NULL, probs = pp_mtx[i, ], method = "Characteristic", log = T)
  }
)

PoissonBinomial_rec_time <- system.time(
  for (i in 1:nrow(pp_mtx)) {
    test_PoissonBinomial_rec[i] <- PoissonBinomial::dpbinom(250, pp_mtx[i, ], method = "Recursive", log = T)
  }
)

# NIMBLE_uc_time <- system.time(
#   for (i in 1:nrow(pp_mtx)) {
#     test_NIMBLE_uc[i] <- dPoisBinom(250, pp_mtx[i, ], log = TRUE)
#   }
# )

NIMBLE_c_time <- system.time(
  for (i in 1:nrow(pp_mtx)) {
    test_NIMBLE_c[i] <- cPoisBinom(250, pp_mtx[i, ], log = TRUE)
  }
)

max(abs(test_NIMBLE_c - test_poisbinom))
max(abs(test_NIMBLE_c - test_NIMBLE_uc))







