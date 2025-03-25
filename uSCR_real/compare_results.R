
# (This doesn't work yet--just moving out of the main script for now)




#### Compare estimated densities ####

result_simple <- readRDS("uSCR_real/simple_posterior.RDS")
simple_area <- diff(c(min(X[, 1]) - 100, max(X[, 1]) + 100)) *
  diff(c(min(X[, 2]) - 100, max(X[, 2]) + 100))

density_samples <- data.frame(
  density = c(
    as.numeric(unlist(result_simple$samples[, "n[1]"])) / simple_area,
    as.numeric(unlist(samples[, "density"])),
    as.numeric(unlist(samples_ROVonly[, "density"]))
  ),
  type = c(
    rep("Cam-only uSCR", length(unlist(result_simple$samples[, "n[1]"]))),
    rep("Joint/masked uSCR", length(unlist(samples[, "density"]))),
    rep("ROV-only GLMM", length(unlist(samples_ROVonly[, "density"])))
  )
)

ggplot(density_samples) +
  geom_density(aes(log(density), group = type, col = type))

# What if we interpret lam0 as part of density

density_samples2 <- data.frame(
  density = c(
    as.numeric(unlist(result_simple$samples[, "n[1]"])) / simple_area* as.numeric(unlist(result_simple$samples[, "lam0"])),
    as.numeric(unlist(samples[, "density"])) * as.numeric(unlist(samples[, "lam0"])),
    as.numeric(unlist(samples_ROVonly[, "density"]))
  ),
  type = c(
    rep("Cam-only uSCR", length(unlist(result_simple$samples[, "n[1]"]))),
    rep("Joint/masked uSCR", length(unlist(samples[, "density"]))),
    rep("ROV-only GLMM", length(unlist(samples_ROVonly[, "density"])))
  )
)
ggplot(density_samples2) +
  geom_density(aes(log(density), group = type, col = type)) +
  xlab("log(density x lam0)")




#### Test points ####
pts <- as.matrix(samples[, c("s[1, 2, 1]", "s[1, 2, 2]")]) %>%
  as.data.frame() %>%
  rename(x = `s[1, 2, 1]`, y = `s[1, 2, 2]`) %>%
  mutate(x = x + x_offset, y = y + y_offset) %>%
  vect(geom = c("x", "y"), crs = crs(hb_mask))

