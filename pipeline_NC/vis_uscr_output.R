library(tidyverse)
library(nimble)
library(gridExtra)

results_files <- list.files("pipeline_NC/NC_results/", 
                            pattern = "uSCR_real_.*.RDS",
                            full.names = TRUE)

summary_df <- lapply(results_files, function(x) {
  temp <- readRDS(x)
  summary <- temp$summary %>% 
    mutate(iter = temp$iter,
           prefix = temp$prefix,
           integration_type = temp$integration_type,
           param = rownames(temp$summary))
  rownames(summary) <- NULL
  summary
}) %>% 
  bind_rows()

inttype_colors <- c(
  "Full" = "#1f78b4",
  "Camera_ROV" = "#b2df8a",
  "Camera_Telemetry" = "#33a02c",
  "Camera_only" = "#fb9a99",
  "Camera_only_noCovar" = "#e31a1c",
  "Telem. prior" = "darkgray",
  "Zulian_etal" = "black"
)

#### Multi-plot visualization - model comparison ####

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
  filter(grepl("^log_sigma", param)) %>% 
  bind_rows(tibble(
    integration_type = "Telem. prior", param = "log_sigma",
    mean = 3.764, `2.5%` = 3.764 - 1.96 * 0.689, `97.5%` = 3.764 + 1.96 * 0.689 
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
for (i in 1:length(results_files)) {
  temp <- readRDS(results_files[i])
  
  ESA_list[[i]] <- calc_ESA(result_list = temp, distr = "Binom") %>% 
    mutate(integration_type = temp$integration_type)
}

plot_ESA <- ESA_list %>% 
  bind_rows() %>% 
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

ggsave("plots/USCR_results.jpg", plot = plot_all, width = 10, height = 6)



#### Multi-plot visualization - full model only ####
inttype_colors_subset <- inttype_colors[c("Full", "Zulian_etal", "Telem. prior")]

# Plot posterior distribution of sigma, 
plot_p0 <- summary_df %>% 
  filter(grepl("p0", param)) %>% 
  filter(integration_type == "Full") %>% 
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
  filter(integration_type == "Full") %>% 
  filter(grepl("^log_sigma", param)) %>% 
  bind_rows(tibble(
    integration_type = "Telem. prior", param = "log_sigma",
    mean = 3.764, `2.5%` = 3.764 - 1.96 * 0.689, `97.5%` = 3.764 + 1.96 * 0.689 
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
  filter(integration_type == "Full") %>% 
  ggplot() +
  geom_pointrange(aes(param, mean/3, ymin = `2.5%`/3, ymax = `97.5%`/3, col = integration_type,
                      group = integration_type), position = position_dodge(width = 0.3)) +
  theme_bw() +
  scale_color_manual("Model", values = inttype_colors) +
  coord_flip() + xlab("") + ggtitle("Population size") +
  ylab("95%CI (number of fish)") + theme(legend.position = "None")

plot_beta <- summary_df %>% 
  filter(integration_type == "Full") %>% 
  filter(grepl("^spatial_beta", param)) %>% 
  filter(integration_type != "Camera_only_noCovar") %>% 
  ggplot() +
  geom_pointrange(aes(param, mean, ymin = `2.5%`, ymax = `97.5%`, col = integration_type,
                      group = integration_type), position = position_dodge(width = 0.3)) +
  theme_bw() +
  scale_color_manual("Model", values = inttype_colors) +
  coord_flip() + xlab("") + ggtitle("Covar. effect") +
  ylab("95%CI")


plot_ESA <- ESA_list %>% 
  bind_rows() %>% 
  filter(integration_type == "Full") %>% 
  ggplot() +
  geom_pointrange(aes(current_dir, ESA_q50, ymin = ESA_q025, ymax = ESA_q975, col = integration_type,
                      group = integration_type), position = position_dodge(width = 0.3)) +
  theme_bw() +
  scale_color_manual("Model", values = inttype_colors) +
  coord_flip() + xlab("Current direction") +
  ggtitle("ESA") +
  ylab("95%CI (m^2)") + theme(legend.position = "None")

plot_ESA_compare <- ESA_list %>% 
  bind_rows() %>% 
  filter(integration_type == "Full") %>% 
  bind_rows(
    data.frame(
      integration_type = "Zulian_etal", current_dir = c("Away", "Across", "Towards"),
      ESA_q50  = c(0.047, 0.064, 0.009) * 1e6, 
      ESA_q025 = c(0.006, 0.012, 0.002) * 1e6, 
      ESA_q975 = c(0.332, 0.346, 0.038) * 1e6
    )
  ) %>%
  mutate(integration_type = factor(integration_type, levels = names(inttype_colors_subset))) %>% 
  ggplot() +
  geom_pointrange(aes(current_dir, ESA_q50, ymin = ESA_q025, ymax = ESA_q975, col = integration_type,
                      group = integration_type), position = position_dodge(width = 0.3),
                  show.legend = TRUE) +
  theme_bw() +
  scale_color_manual("Model", values = inttype_colors_subset, drop = FALSE) +
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

ggsave("plots/USCR_results_fullOnly.jpg", plot = plot_all, width = 10, height = 6)

