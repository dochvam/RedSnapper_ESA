library(nimble)
library(tidyverse)

temp <- readRDS("pipeline_NC/NC_results/uSCR_real_Augustine_Binom_Final_100_Full.RDS")

target_cols <- c(
  paste0("s[", 1:3000, ", 1]"),
  paste0("s[", 1:3000, ", 2]"),
  paste0("z[", 1:3000, "]")
)
df <- temp$samples2 %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  select(all_of(target_cols))


vis_one_post <- function(df, target_ID) {
  this_df <- df %>% 
    dplyr::select(all_of(c(
      paste0("s[", target_ID, ", 1]"),
      paste0("s[", target_ID, ", 2]"),
      paste0("z[", target_ID, "]")
    )))
  colnames(this_df) <- c("x", "y", "z")
  
  s_plot <- ggplot(this_df) + 
    geom_point(aes(x, y, col = as.factor(z))) +
    theme_minimal() +
    scale_color_manual(values = c("pink", "black")) +
    ggtitle(paste0("Individual index ", target_ID))
  
  s_plot
}

vis_one_post(df, 100)
vis_one_post(df, 101)
vis_one_post(df, 1000)
vis_one_post(df, 1002)
vis_one_post(df, 1902)


#### Visualize all of them together #### 

s1_cols <- paste0("s[", 1:600, ", 1]")
s2_cols <- paste0("s[", 1:600, ", 2]")
z_cols <- paste0("z[", 1:600, "]")

s1_long <- df[, s1_cols] %>% 
  pivot_longer(cols = all_of(s1_cols)) %>% 
  rename(s1 = value)
s2_long <- df[, s2_cols] %>% 
  pivot_longer(cols = all_of(s2_cols)) %>% 
  rename(s2 = value)
z_long <- df[, z_cols] %>% 
  pivot_longer(cols = all_of(z_cols)) %>% 
  rename(z = value)

all_obs_df <- bind_cols(
  s1 = s1_long$s1, s2 = s2_long$s2, z = z_long$z
)

all_obs_df %>% 
  filter(z == 0) %>% 
  ggplot() +
  geom_bin2d(aes(s1, s2))

#### Visualize the posterior distribution of covar surface effect #### 
spatial_beta <- mean(temp$samples[, "spatial_beta"])
phi <- exp(log(vps_mtx) * spatial_beta)
image(phi / sum(phi))
