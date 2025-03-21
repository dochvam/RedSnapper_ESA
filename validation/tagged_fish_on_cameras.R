library(readxl)
library(tidyverse)
library(terra)

all_stations <-  read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/cameratraps/RS_2023_full_reads_all_three_cameratrap_dates.xlsx", sheet = "StationData") %>% 
  distinct(Station_ID, Date, 
           Longitude = Start_Longitude, 
           Latitude = Start_Latitude) %>% 
  mutate(Longitude = as.numeric(Longitude))

tagged_fish_cameras <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/cameratraps/RS_2023_full_reads_all_three_cameratrap_dates.xlsx") %>% 
  filter(grepl("Vemco", comments)) %>% 
  mutate(fish_tag = paste0(125, parse_number(comments))) %>% 
  distinct(fish_tag, Station_ID, Date = date) %>% 
  left_join(
   all_stations
  )

for (i in 1:nrow(tagged_fish_cameras)) {
  this_fish <- tagged_fish_cameras$fish_tag[i]
  this_ID <- tagged_fish_cameras$Station_ID[i]
    
  this_fish_vps <- read_csv(paste0(
    "Data/SnapperMvmtAbundanceStudy/VPS_Data/VPS-ChickenRock-01-Results-20240202/results/animal/", 
    this_fish, ".csv")
  ) %>% 
    filter(as_date(Time) == tagged_fish_cameras$Date[i])
  
  today_stations <- all_stations %>% 
    filter(Date == tagged_fish_cameras$Date[i])
  
  this_plot <- ggplot(mapping = aes(Longitude, Latitude)) +
    geom_point(data = this_fish_vps, col = "blue") + 
    geom_point(data = today_stations, 
               aes(col = Station_ID == this_ID),
               show.legend = F,
               pch = 18, size = 2) +
    theme_minimal() + 
    ggtitle(paste0("Fish ", this_fish)) +
    scale_color_manual("", values = c("black", "red"))
  ggsave(paste0("plots/fish_on_cameras/", this_fish, ".jpg"), this_plot,
         width = 4, height = 3)
}
