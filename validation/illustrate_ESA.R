library(tidyverse)


resoln <- 5
max_m <- 500
cell_area <- resoln^2
grid <- expand.grid(x = seq(-max_m, max_m, by = resoln),
                    y = seq(-max_m, max_m, by = resoln)) %>% 
  as.data.frame() %>% 
  mutate(dist = sqrt(x^2 + y^2))



df <- expand.grid(lam0 = c(0.5, 0.9), sigma = c(50, 100))

outlist <- list()
for (i in 1:nrow(df)) {
  outlist[[i]] <- grid %>% 
    mutate(lam0 = df$lam0[i], sigma = df$sigma[i])
  outlist[[i]]$lambda <- df$lam0[i] * exp(-grid$dist^2 / (2 * df$sigma[i]^2))
  outlist[[i]]$prob_one_det <- 1-exp(-outlist[[i]]$lambda)
}

bind_rows(outlist) %>% 
  ggplot() +
  geom_tile(aes(x, y, fill = prob_one_det)) +
  facet_grid(paste0("sigma = ", sigma)~paste0("lambda0 = ", lam0)) +
  ggtitle("ESA illustration for Poisson response model") +
  theme_bw() +
  scale_fill_viridis_c()


