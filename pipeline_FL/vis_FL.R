library(tidyverse)
library(nimble)
library(gridExtra)
library(coda)
library(MCMCvis)
library(gridExtra)

results_files_full_HBmap <- list.files("pipeline_FL/FL_results/", 
                                 pattern = "uSCR_real_Augustine_Binom_*.*Full_HB_map.RDS",
                                 full.names = TRUE)

samples_full_HBmap <- lapply(results_files_full_HBmap, function(x) readRDS(x)$samples) %>% 
  as.mcmc.list()

summary_full_HBmap <- MCMCsummary(samples_full_HBmap) %>% 
  mutate(integration_type = "Full", spatial_type = "HB_map")
summary_full_HBmap$param <- rownames(summary_full_HBmap)
rownames(summary_full_HBmap) <- NULL

results_files_full_VPSmap <- list.files("pipeline_FL/FL_results/", 
                                 pattern = "uSCR_real_Augustine_Binom_*.*Full_VPS_map.RDS",
                                 full.names = TRUE)

samples_full_VPSmap <- lapply(results_files_full_VPSmap, function(x) readRDS(x)$samples) %>% 
  as.mcmc.list()

summary_full_VPSmap <- MCMCsummary(samples_full_VPSmap) %>% 
  mutate(integration_type = "Full", spatial_type = "VPS_map")
summary_full_VPSmap$param <- rownames(summary_full_VPSmap)
rownames(summary_full_VPSmap) <- NULL


results_files_camera_telemetry_HBmap <- list.files("pipeline_FL/FL_results/", 
                                 pattern = "uSCR_real_Augustine_Binom_*.*Camera_Telemetry_HB_map.RDS",
                                 full.names = TRUE)

samples_camera_telemetry_HBmap <- lapply(results_files_camera_telemetry_HBmap, function(x) readRDS(x)$samples) %>% 
  as.mcmc.list()

summary_camera_telemetry_HBmap <- MCMCsummary(samples_camera_telemetry_HBmap) %>% 
  mutate(integration_type = "Camera_Telemetry", spatial_type = "HB_map")
summary_camera_telemetry_HBmap$param <- rownames(summary_camera_telemetry_HBmap)
rownames(summary_camera_telemetry_HBmap) <- NULL

results_files_camera_telemetry_VPSmap <- list.files("pipeline_FL/FL_results/", 
                                 pattern = "uSCR_real_Augustine_Binom_*.*Camera_Telemetry_VPS_map.RDS",
                                 full.names = TRUE)

samples_camera_telemetry_VPSmap <- lapply(results_files_camera_telemetry_VPSmap, function(x) readRDS(x)$samples) %>% 
  as.mcmc.list()

summary_camera_telemetry_VPSmap <- MCMCsummary(samples_camera_telemetry_VPSmap) %>% 
  mutate(integration_type = "Camera_Telemetry", spatial_type = "VPS_map")
summary_camera_telemetry_VPSmap$param <- rownames(summary_camera_telemetry_VPSmap)
rownames(summary_camera_telemetry_VPSmap) <- NULL


# results_files_cam <- list.files("pipeline_NC/NC_results/", 
#                                 pattern = "uSCR_real_Augustine_Binom_20minSigma_ScalePrior*.*Camera_only.RDS",
#                                 full.names = TRUE)
# 
# samples_cam <- lapply(results_files_cam, function(x) readRDS(x)$samples) %>% 
#   as.mcmc.list()
# 
# summary_cam <- MCMCsummary(samples_cam) %>% 
#   mutate(integration_type = "Camera_only")
# summary_cam$param <- rownames(summary_cam)
# rownames(summary_cam) <- NULL

summary_df <- bind_rows(summary_full_HBmap, summary_full_VPSmap,
                        summary_camera_telemetry_HBmap, summary_camera_telemetry_VPSmap) %>% 
  mutate(integration_type = paste0(integration_type, "_", spatial_type))


#### Multi-plot visualization - model comparison ####
# inttype_colors <- c(
#   "Full" = "#1f78b4",
#   "Camera_only" = "#e31a1c",
#   "Telem. prior" = "darkgray",
#   "Zulian_etal" = "black"
# )
inttype_colors <- c(
  "Full_HB_map" = "#1f78b4",
  "Full_VPS_map" = "darkgreen",
  "Camera_Telemetry_HB_map" = "#e31a1c",
  "Camera_Telemetry_VPS_map" = "darkred",
  "Telem. prior" = "darkgray",
  "Zulian_etal" = "black"
)
log_sigma_estimate <- read_csv("pipeline_FL/FL_results/log_sigma_estimate_FL.csv")
log_sigma_mean <- 1.93
log_sigma_sd <- 1.02


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
ESA_list[[1]] <- calc_ESA(result_list = list(samples = samples_full_HBmap), distr = "Binom") %>% 
  mutate(integration_type = "Full_HB_map")
ESA_list[[2]] <- calc_ESA(result_list = list(samples = samples_full_VPSmap), distr = "Binom") %>% 
  mutate(integration_type = "Full_VPS_map")
ESA_list[[3]] <- calc_ESA(result_list = list(samples = samples_camera_telemetry_HBmap), distr = "Binom") %>% 
  mutate(integration_type = "Camera_Telemetry_HB_map")
ESA_list[[4]] <- calc_ESA(result_list = list(samples = samples_camera_telemetry_VPSmap), distr = "Binom") %>% 
  mutate(integration_type = "Camera_Telemetry_VPS_map")

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

ggsave("plots/FL_USCR_20min_results.jpg", plot = plot_all, width = 10, height = 6)


