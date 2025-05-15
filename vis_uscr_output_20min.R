library(tidyverse)
library(nimble)
library(gridExtra)
library(coda)
library(MCMCvis)
library(gridExtra)

results_files_full <- list.files("pipeline_NC/NC_results/", 
                            pattern = "uSCR_real_Augustine_Binom_20minSigma_ScalePrior*.*Full.RDS",
                            full.names = TRUE)

samples_full <- lapply(results_files_full, function(x) readRDS(x)$samples) %>% 
  as.mcmc.list()

summary_full <- MCMCsummary(samples_full) %>% 
  mutate(integration_type = "Full")
summary_full$param <- rownames(summary_full)
rownames(summary_full) <- NULL


results_files_cam <- list.files("pipeline_NC/NC_results/", 
                            pattern = "uSCR_real_Augustine_Binom_20minSigma_ScalePrior*.*Camera_only.RDS",
                            full.names = TRUE)

samples_cam <- lapply(results_files_cam, function(x) readRDS(x)$samples) %>% 
  as.mcmc.list()

summary_cam <- MCMCsummary(samples_cam) %>% 
  mutate(integration_type = "Camera_only")
summary_cam$param <- rownames(summary_cam)
rownames(summary_cam) <- NULL

summary_df <- bind_rows(summary_full, summary_cam)


#### Multi-plot visualization - model comparison ####
inttype_colors <- c(
  "Full" = "#1f78b4",
  "Camera_only" = "#e31a1c",
  "Telem. prior" = "darkgray",
  "Zulian_etal" = "black"
)
# log_sigma_estimate <- read_csv("pipeline_NC/NC_results/log_sigma_estimate_NC.csv")
log_sigma_mean <- 2.59
log_sigma_sd <- 0.803


# Plot posterior distribution of sigma, 
plot_p0 <- summary_df %>% 
  filter(grepl("p0", param)) %>% 
  mutate(dir_num = parse_number(substr(param,4,4)),
         dir = c("Towards", "Away", "Across")[dir_num]) %>% 
  ggplot() +
  geom_pointrange(aes(dir, mean, ymin = `2.5%`, ymax = `97.5%`, col = integration_type,
                      group = integration_type), position = position_dodge(width = 0.3)) +
  theme_bw() +
  scale_color_manual("Model", values = inttype_colors) +
  coord_flip() + xlab("") + ggtitle("Baseline detection prob.") +
  ylab("95%CI") + theme(legend.position = "None")

plot_sigma <- summary_df %>% 
  filter(grepl("^sigma", param)) %>% 
  bind_rows(tibble(
    integration_type = "Telem. prior", param = "sigma",
    mean = exp(log_sigma_mean), 
    `2.5%` = exp(log_sigma_mean - 1.96 * log_sigma_sd),
    `97.5%` = exp(log_sigma_mean + 1.96 * log_sigma_sd)
  )) %>% 
  ggplot() +
  geom_pointrange(aes(param, mean, ymin = `2.5%`, ymax = `97.5%`, col = integration_type,
                      group = integration_type), position = position_dodge(width = 0.3)) +
  theme_bw() +
  scale_color_manual("Model", values = inttype_colors) +
  coord_flip() + xlab("") + ggtitle("Log-scale displacement param.") +
  ylab("95%CI (m)") + theme(legend.position = "None")


plot_N <- summary_df %>% 
  filter(grepl("^N", param)) %>% 
  ggplot() +
  geom_pointrange(aes(param, mean/3, ymin = `2.5%`/3, ymax = `97.5%`/3, col = integration_type,
                      group = integration_type), position = position_dodge(width = 0.3)) +
  theme_bw() +
  scale_color_manual("Model", values = inttype_colors) +
  coord_flip() + xlab("") + ggtitle("Population size") +
  ylab("95%CI (number of fish)") + theme(legend.position = "None")

plot_beta <- summary_df %>% 
  filter(grepl("^spatial_beta", param)) %>% 
  filter(integration_type != "Camera_only_noCovar") %>% 
  ggplot() +
  geom_pointrange(aes(param, mean, ymin = `2.5%`, ymax = `97.5%`, col = integration_type,
                      group = integration_type), position = position_dodge(width = 0.3)) +
  theme_bw() +
  scale_color_manual("Model", values = inttype_colors) +
  coord_flip() + xlab("") + ggtitle("Covar. effect") +
  ylab("95%CI")



source("pipeline_NC/ESA_helper.R")

ESA_list <- list()
ESA_list[[1]] <- calc_ESA(result_list = list(samples = samples_cam), distr = "Binom") %>% 
  mutate(integration_type = "Camera_only")
ESA_list[[2]] <- calc_ESA(result_list = list(samples = samples_full), distr = "Binom") %>% 
  mutate(integration_type = "Full")

plot_ESA <- ESA_list %>% 
  bind_rows() %>% 
  filter(integration_type != "Camera_only_noCovar") %>% 
  ggplot() +
  geom_pointrange(aes(current_dir, ESA_q50, ymin = ESA_q025, ymax = ESA_q975, col = integration_type,
                      group = integration_type), position = position_dodge(width = 0.3)) +
  theme_bw() +
  scale_color_manual("Model", values = inttype_colors) +
  coord_flip() + xlab("") + 
  ggtitle("ESA") +
  ylab("95%CI (m^2)") + theme(legend.position = "None")

plot_ESA_compare <- ESA_list %>% 
  bind_rows() %>% 
  filter(integration_type != "Camera_only_noCovar") %>% 
  bind_rows(
    data.frame(
      integration_type = "Zulian_etal", current_dir = c("Away", "Across", "Towards"),
      ESA_q50  = c(0.047, 0.064, 0.009) * 1e6, 
      ESA_q025 = c(0.006, 0.012, 0.002) * 1e6, 
      ESA_q975 = c(0.332, 0.346, 0.038) * 1e6
    )
  ) %>%
  mutate(integration_type = factor(integration_type, levels = names(inttype_colors))) %>% 
  ggplot() +
  geom_pointrange(aes(current_dir, ESA_q50, ymin = ESA_q025, ymax = ESA_q975, col = integration_type,
                      group = integration_type), position = position_dodge(width = 0.3),
                  show.legend = TRUE) +
  theme_bw() +
  scale_color_manual("Model", values = inttype_colors, drop = FALSE) +
  coord_flip() + xlab("") + 
  ggtitle("ESA (incl. previous estimate)") +
  ylab("95%CI (m^2)")

modeltype_legend <- get_legend(plot_ESA_compare)
plot_ESA_compare <- plot_ESA_compare + theme(legend.position = "None")

layout <- matrix(
  c(1, 4,
    2, 4,
    3, 5,
    3, 5),
  ncol = 2, byrow = TRUE
)

plot_all <- arrangeGrob(plot_N, plot_sigma, plot_p0, plot_ESA, 
                        plot_ESA_compare, layout_matrix = layout,
                        right = modeltype_legend)

ggsave("plots/USCR_20min_results.jpg", plot = plot_all, width = 10, height = 6)


