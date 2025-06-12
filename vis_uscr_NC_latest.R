library(tidyverse)
library(nimble)
library(gridExtra)
library(coda)
library(MCMCvis)
library(gridExtra)

source("pipeline_NC/ESA_helper.R")

results_files_all <- list.files("pipeline_NC/NC_results/", 
                                pattern = "uSCR_real_Augustine_Binom_FinalChoices*.*RDS",
                                full.names = TRUE)
results_files_all <- results_files_all[grepl("_Offset", results_files_all)]

combos <- expand.grid(
  # integration_type = c("Full_Offset", "Full_NoOffset", "Camera_Only"),
  integration_type = c("Full_Offset"),
  # sigma_type = c("Constant", "Prior", "Trunc_prior")
  sigma_type = c("Prior_Mean", "Prior_Variability")
)


samples_list_list <- list()
summary_list <- list()
for (i in 1:nrow(combos)) {
  this_resfiles <- results_files_all[grepl(paste0(combos$integration_type[i], "_", combos$sigma_type[i]),
                                           results_files_all)]
  
  samples_list_list[[i]] <- lapply(this_resfiles, function(x) readRDS(x)$samples) %>% 
    as.mcmc.list()
  summary_list[[i]] <- MCMCsummary(samples_list_list[[i]]) %>%
    mutate(integration_type = combos$integration_type[i],
           sigma_type = combos$sigma_type[i])
  summary_list[[i]]$param <- rownames(summary_list[[i]])
  rownames(summary_list[[i]]) <- NULL
  
}



summary_df <- bind_rows(summary_list) %>% 
  mutate(type = paste0(integration_type, "_", sigma_type))


#### Multi-plot visualization - model comparison ####
type_colors <- c(
  viridisLite::plasma(length(unique(combos$sigma_type)), end = 0.9), "darkgray", "black"
)
names(type_colors) <- c(
  unique(as.character(combos$sigma_type)),
  "Prior",
  "Zulian et al. ESA"
)

capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

log_sigma_estimate <- read_csv("pipeline_NC/NC_results/log_sigma_estimate_NC.csv") %>% 
  mutate(sigma_type = capwords(type), integration_type = "Prior")

# Plot posterior distribution of sigma, 
plot_p0 <- summary_df %>% 
  filter(grepl("p0", param)) %>% 
  mutate(dir_num = parse_number(substr(param,4,4)),
         dir = c("Towards", "Away", "Across")[dir_num]) %>% 
  ggplot() +
  geom_pointrange(aes(dir, mean, ymin = `2.5%`, ymax = `97.5%`, col = sigma_type,
                      group = type, shape = sigma_type), position = position_dodge(width = 0.4)) +
  theme_bw() +
  scale_color_manual("Model", values = type_colors) +
  coord_flip() + xlab("") + ggtitle("Baseline detection prob.") +
  ylab("95%CI") + theme(legend.position = "None") +
  scale_shape_manual("Model", values = c(1, 17, 19))


log_sigma_estimate_toplot <- log_sigma_estimate %>% 
  filter(type == "mean")

attraction_dist_sigma <- 0
prior_df <- tibble(
  integration_type = "Prior", param = "sigma",
  sigma_type = c(log_sigma_estimate_toplot$sigma_type),
  type = paste0("Prior_", log_sigma_estimate_toplot$sigma_type),
  mean = sqrt(exp(log_sigma_estimate_toplot$mean)^2 + attraction_dist_sigma^2),
  `2.5%` = sqrt(exp(log_sigma_estimate_toplot$mean - 1.96 * log_sigma_estimate_toplot$sd)^2 + attraction_dist_sigma^2),
  `97.5%` = sqrt(exp(log_sigma_estimate_toplot$mean + 1.96 * log_sigma_estimate_toplot$sd)^2 + attraction_dist_sigma^2)
)

plot_sigma <- summary_df %>% 
  filter(grepl("^sigma", param)) %>% 
  bind_rows(prior_df) %>% 
  ggplot() +
  geom_pointrange(aes(param, mean, ymin = `2.5%`, ymax = `97.5%`, col = sigma_type,
                      group = type, shape = sigma_type), position = position_dodge(width = 0.9)) +
  theme_bw() +
  scale_color_manual("Model", values = type_colors) +
  coord_flip() + xlab("") + ggtitle("Displacement param. (m)") +
  ylab("95%CI (m)") +
  scale_shape_manual("Model", values = c(1, 1, 1, 18, 19))
modeltype_legend <- get_legend(plot_sigma)
plot_sigma <- plot_sigma + theme(legend.position = "None")


plot_N <- summary_df %>% 
  filter(grepl("^N", param)) %>% 
  ggplot() +
  geom_pointrange(aes(param, mean/3, ymin = `2.5%`/3, ymax = `97.5%`/3, col = sigma_type,
                      group = type, shape = sigma_type), position = position_dodge(width = 0.4)) +
  theme_bw() +
  scale_color_manual("Model", values = type_colors) +
  coord_flip() + xlab("") + ggtitle("Population size") +
  ylab("95%CI (number of fish)") + theme(legend.position = "None") +
  scale_shape_manual("Model", values = c(1, 18, 19))

plot_beta <- summary_df %>% 
  filter(grepl("^spatial_beta", param)) %>% 
  filter(integration_type != "Camera_only_noCovar") %>% 
  ggplot() +
  geom_pointrange(aes(param, mean, ymin = `2.5%`, ymax = `97.5%`, col = sigma_type,
                      group = type, shape = sigma_type), position = position_dodge(width = 0.4)) +
  theme_bw() +
  scale_color_manual("Model", values = type_colors) +
  coord_flip() + xlab("") + ggtitle("Covar. effect") +
  ylab("95%CI") +
  scale_shape_manual("Model", values = c(1, 18, 19))



esa_types <- c("20-min", 
               "NB attraction dist",
               "NB attraction dist plus 20-min",
               "as modeled")

ESA_list <- list()
ct <- 0
for (i in 1:nrow(combos)) {
  for (j in 1:length(esa_types)) {
    ct <- ct + 1
    ESA_list[[ct]] <- calc_ESA(result_list = list(samples = samples_list_list[[i]]), distr = "Binom",
                               log_sigma_estimate = log_sigma_estimate,
                               esa_sigma_type = esa_types[j]) %>% 
      mutate(integration_type = summary_list[[i]]$integration_type[1], 
             sigma_type = summary_list[[i]]$sigma_type[1],
             esa_type = esa_types[j],
             type = paste0(integration_type, "_", sigma_type))
  }
}

ESA_list %>% 
  bind_rows() %>% 
  write_csv("Chicken_Rock_ESA_table_0603.csv")

plot_ESA <- ESA_list %>% 
  bind_rows() %>% 
  # filter(sigma_type == "Prior_Variability") %>% 
  ggplot() +
  geom_pointrange(aes(current_dir, ESA_q50, ymin = ESA_q025, ymax = ESA_q975, 
                      col = esa_type,
                      group = esa_type), position = position_dodge(width = 0.4)) +
  theme_bw() +
  scale_color_viridis_d("Sigma type", option = "turbo") +
  coord_flip() + xlab("") + 
  ggtitle("ESA") +
  ylab("95%CI (m^2)") +
  facet_wrap(~sigma_type, ncol = 1) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2))

plot_ESA_compare <- ESA_list %>% 
  bind_rows() %>% 
  bind_rows(
    data.frame(
      type = "Zulian et al. ESA", 
      sigma_type = "Zulian et al. ESA",
      current_dir = c("Away", "Across", "Towards"),
      ESA_q50  = c(0.047, 0.064, 0.009) * 1e6, 
      ESA_q025 = c(0.006, 0.012, 0.002) * 1e6, 
      ESA_q975 = c(0.332, 0.346, 0.038) * 1e6
    )
  ) %>%
  # mutate(integration_type = factor(integration_type, levels = names(type_colors))) %>% 
  ggplot() +
  geom_pointrange(aes(current_dir, ESA_q50, ymin = ESA_q025, ymax = ESA_q975, col = sigma_type,
                      group = type, shape = sigma_type), position = position_dodge(width = 0.4),
                  show.legend = TRUE) +
  theme_bw() +
  scale_color_manual("Model", values = type_colors, drop = FALSE) +
  coord_flip() + xlab("") + 
  ggtitle("ESA (incl. previous estimate)") +
  ylab("95%CI (m^2)") +
  scale_shape_manual("Model", values = c(1, 18, 19, 19))


layout <- matrix(
  c(1, 4,
    2, 4,
    3, 4,
    3, 4),
  ncol = 2, byrow = TRUE
)

plot_all <- arrangeGrob(plot_N, plot_sigma, plot_p0, plot_ESA, 
                        layout_matrix = layout,
                        left = modeltype_legend)

ggsave("plots/USCR_20min_results.jpg", plot = plot_all, width = 11, height = 6)


