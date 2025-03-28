library(tidyverse)
library(readxl)

n <- 32

camera_dat <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/cameratraps/RS_2023_full_reads_all_three_cameratrap_dates.xlsx") %>% 
  mutate(time = as.numeric(difftime(time, as_datetime("1899-12-31 00:00:00"), units = "sec"))) %>% 
  group_by(Station_ID, date) %>% 
  summarize(total = max(total),
            tagged = max(tagged)) %>% 
  filter(!is.na(Station_ID))

k_tagged <- sum(camera_dat$tagged)
k_total <- sum(camera_dat$total)

# Chapman estimator
(n + 1) * (k_total + 1) / (k_tagged + 1)

# LP estimator
(n) * (k_total) / (k_tagged)
