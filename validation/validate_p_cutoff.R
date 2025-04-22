
cutoff_df <- data.frame(
  p_cutoff = 10^-(2:11),
  n_excluded = NA,
  ll_error = NA,
  time = NA
)

cmod$p_cutoff <- 1e-100
full_ll <- cmod$calculate()

for (i in 1:nrow(cutoff_df)) {
  cutoff_df$n_excluded[i] <- sum(cmod$detprob < cutoff_df$p_cutoff[i], na.rm = T)
  cmod$p_cutoff <- cutoff_df$p_cutoff[i]
  cutoff_df$time[i] <- system.time(cutoff_df$ll_error[i] <- cmod$calculate() - full_ll)
}


p1 <- ggplot(cutoff_df) +
  geom_line(aes(as.factor(p_cutoff), ll_error, group = 1)) + ylab("Error in log-likelihood calc.") +
  xlab("P cutoff") + theme_minimal()

p2 <- ggplot(cutoff_df) +
  geom_line(aes(as.factor(p_cutoff), time, group = 1)) + ylab("Computation time") +
  xlab("P cutoff") + theme_minimal()

p3 <- ggplot(cutoff_df) +
  geom_line(aes(as.factor(p_cutoff), n_excluded, group = 1)) + ylab("Num. fish-cam pairs excluded") +
  xlab("P cutoff") + theme_minimal()


p_all <- gridExtra::arrangeGrob(p1, p2, p3, nrow = 1)
ggsave(plot = p_all, "plots/validate_pcutoff.jpg", width = 12, height = 3)


