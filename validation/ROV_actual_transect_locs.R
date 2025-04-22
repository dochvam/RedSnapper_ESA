library(readxl)
library(tidyverse)
library(terra)
# install.packages("measurements") 
library(measurements)


temp <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/ROV/ROV lat and long locations for transects in 2023.xlsx", col_names = F)
colnames(temp) <- paste0("X", 1:5)

temp$X5[159] <- "34 37.357N 76 08.686W"
temp$X5[63] <-  "34 37.437N 76 08.686W"

for (i in c(36:39, 51:60, 126:134)) {
  for (col in paste0("X", 3:5)) {
    temp[i, col] <- gsub(" 08", " 09", temp[i, col])
  }
}



reprocess <- distinct(temp, X1, X2) %>% 
  rename(Site_ID = X1, Date_Numeric = X2) %>% 
  filter(Site_ID != "Site:", !is.na(Site_ID)) %>% 
  mutate(Date = as_date(ifelse(Date_Numeric == min(Date_Numeric), "2023-08-21", "2023-08-22"))) %>% 
  mutate(
    drop1_ll = NA,
    start1_ll = NA,
    end1_ll = NA,
    drop2_ll = NA,
    start2_ll = NA,
    end2_ll = NA
  )

for (i in 1:nrow(reprocess)) {
  this_start_ind <- which(temp$X1 == reprocess$Site_ID[i] & temp$X2 == as.numeric(reprocess$Date_Numeric[i]))
  reprocess$drop1_ll[i]  <- temp$X4[this_start_ind]
  reprocess$start1_ll[i] <- temp$X5[this_start_ind - 1]
  reprocess$end1_ll[i]   <- temp$X5[this_start_ind]
  reprocess$drop2_ll[i]  <- temp$X4[this_start_ind + 2]
  reprocess$start2_ll[i] <- temp$X5[this_start_ind + 1]
  reprocess$end2_ll[i]   <- temp$X5[this_start_ind + 2]
  # reprocess$drop1_ll[i]  <- temp$X5[this_start_ind - 1]
  # reprocess$start1_ll[i] <- temp$X4[this_start_ind]
  # reprocess$end1_ll[i]   <- temp$X5[this_start_ind]
  # reprocess$drop2_ll[i]  <- temp$X5[this_start_ind + 1]
  # reprocess$start2_ll[i] <- temp$X4[this_start_ind + 2]
  # reprocess$end2_ll[i]   <- temp$X5[this_start_ind + 2]
}
    

ROV_transect_locs <- reprocess %>% 
  pivot_longer(cols = c("drop1_ll", "start1_ll", "end1_ll", "drop2_ll", "start2_ll", "end2_ll")) %>% 
  mutate(transect = parse_number(name), point = substr(name, 1, nchar(name) - 4)) %>% 
  filter(value != "none", value != "NA") %>% 
  # filter(point %in% c("start", "end")) %>% 
  mutate(latitude = as.numeric(conv_unit(substr(value, 1, 9), from = "deg_dec_min", to = "dec_deg")),
         longitude = -1 * as.numeric(conv_unit(substr(value, 12, 20), from = "deg_dec_min", to = "dec_deg"))) %>% 
  select(Site_ID, Date, point, transect, longitude, latitude) %>% 
  group_by(Site_ID, Date, transect) %>% 
  mutate(point = factor(point, levels = c("start", "drop", "end"))) %>% 
  arrange(Date, Site_ID, transect, point)
  # mutate(ll_dist = sqrt(diff(longitude)^2 + diff(latitude)^2))



write_csv(ROV_transect_locs, "Data/SnapperMvmtAbundanceStudy/CountData/ROV/ROV_latlong_real_2023_processed.csv")



ggplot(ROV_transect_locs) +
  geom_path(aes(longitude, latitude, group = paste0(Site_ID, transect, Date))) +
  # geom_point(aes(longitude, latitude, col = point)) +
  facet_wrap(~Date, nrow = 2) +
  theme_bw() + coord_fixed() # +
  # geom_label(aes(longitude, latitude, label = Site_ID), size = 2)


