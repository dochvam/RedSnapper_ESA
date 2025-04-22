library(tidyverse)

esa_files <- list.files("ESA_estimates/", pattern = "ESA.csv",
                        full.names = TRUE)

esa_results <- lapply(esa_files, read_csv) %>% 
  bind_rows() %>% 
  mutate(current_dir = ifelse(is.na(current_dir), dir, current_dir)) %>% 
  select(-dir) %>% 
  rename(type = prefix)

esa_results$type[grepl("VPSsurface", esa_results$type)] <- "uSCR_VPS_surface"
esa_results$type[grepl("VPSasCovar_uSCR_sansROV", esa_results$type)] <- "uSCR_VPS_asCovar_sansROV"
esa_results$type[grepl("VPSasCovar_uSCR_wROV", esa_results$type)] <- "uSCR_VPS_asCovar_wROV"
esa_results$type[grepl("VPSasCovar_uninformative", esa_results$type)] <- "uSCR_VPS_asCovar_uninformativePrior"
esa_results$type[grepl("masked", esa_results$type)] <- "uSCR_hab_mask"
esa_results$name <- paste0(esa_results$type, "_", esa_results$distr)

colors <- c(
  # "#68956E",
  "#226F54",
  "#87C38F",
  "#DA2C38",
  "#711218",
  "#43291F",
  "black"
  # "#8C746B",
  # "#cfa495"
)


# esa_results %>%
#   filter(!grepl("Count", type)) %>%
#   bind_rows(data.frame(
#     name = "Zulian_etal",
#     ESA_q50 = 0.047 * 1e6, ESA_q025 = 0.006 * 1e6, ESA_q975 = 0.332 * 1e6
#   )) %>%
#   ggplot() +
#   geom_pointrange(aes(current_dir, ESA_q50,
#                       ymin = ESA_q025, ymax = ESA_q975,
#                       col = name, shape = name), position = position_dodge(width = 0.3)) +
#   coord_flip() +
#   theme_minimal() +
#   scale_color_manual("Model", values = colors) +
#   scale_shape_manual("Model", values = c(1, 19, 24, 1, 19, 1, 19, 24, 25)) +
#   scale_y_continuous(trans = "log", breaks = 2^(7:13)) +
#   ylab("ESA 95%CI (in m^2; note log-scale axis)") + xlab("")


ESA_compare <- esa_results %>%
  filter(!grepl("Count", type)) %>%
  bind_rows(data.frame(
    name = "Zulian_etal", current_dir = c("Away", "Perpendicular", "Towards"),
    ESA_q50  = c(0.047, 0.064, 0.009) * 1e6, 
    ESA_q025 = c(0.006, 0.012, 0.002) * 1e6, 
    ESA_q975 = c(0.332, 0.346, 0.038) * 1e6
  )) %>%
  ggplot() +
  geom_pointrange(aes(name, ESA_q50, 
                      ymin = ESA_q025, ymax = ESA_q975, 
                      col = current_dir, shape = current_dir), position = position_dodge(width = 0.3)) +
  coord_flip() +
  scale_color_viridis_d(end = 0.8) +
  theme_minimal() +
  scale_y_continuous(trans = "log", breaks = 4^(3:9)) +
  ylab("ESA 95%CI (in m^2; log-scale x-axis)") + xlab("")

esa_results %>%
  filter(!grepl("Count", type)) %>%
  bind_rows(data.frame(
    name = "Zulian_etal", current_dir = c("Away", "Perpendicular", "Towards"),
    ESA_q50  = c(0.047, 0.064, 0.009) * 1e6, 
    ESA_q025 = c(0.006, 0.012, 0.002) * 1e6, 
    ESA_q975 = c(0.332, 0.346, 0.038) * 1e6
  )) %>%
  ggplot() +
  geom_pointrange(aes(name, ESA_q50, 
                      ymin = ESA_q025, ymax = ESA_q975, 
                      col = current_dir, shape = current_dir), position = position_dodge(width = 0.3)) +
  coord_flip() +
  scale_color_viridis_d(end = 0.8) +
  theme_minimal() +
  # scale_y_continuous(trans = "log", breaks = 4^(3:9)) +
  ylab("ESA 95%CI (in m^2; real scale)") + xlab("")

ggsave(plot = ESA_compare, "plots/ESA_comparison.jpg",
       width = 6, height = 4, dpi = 300)


