library(tidyverse)
library(terra)


#### Load the camera data ####
M <- 1000

# Get the real camera locations and deployment times for use in the simulation

camera_dat <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/cameratraps/RS_2023_full_reads_all_three_cameratrap_dates.xlsx") %>% 
  mutate(time = as.numeric(difftime(time, as_datetime("1899-12-31 00:00:00"), units = "sec"))) %>% 
  filter(camera == "A", time >= 60 * 10) # Filter out descent/retrieval frames

camera_counts <- camera_dat %>%
  group_by(Station_ID, date) %>% 
  summarize(count = max(total))

camera_locs <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/cameratraps/RS_2023_full_reads_all_three_cameratrap_dates.xlsx", sheet = "StationData") %>% 
  filter(`Camera (A or C)` == "A") %>% 
  dplyr::select(Station_ID, Date, Time = Start_Time_GMT,
                Latitude = Start_Latitude, Longitude = Start_Longitude) %>% 
  left_join(camera_counts, by = c("Station_ID", "Date" = "date"))


camera_locs$Longitude <- as.numeric(camera_locs$Longitude)

camera_pts <- vect(camera_locs, geom = c("Longitude", "Latitude"),
                   crs = "+proj=longlat") %>% 
  project("ESRI:102003") %>% 
  as.data.frame(geom = "XY")

X <- camera_pts[, c("x", "y")]
X[, 1] <- X[, 1] - mean(X[, 1])
X[, 2] <- X[, 2] - mean(X[, 2])
xlim <- c(min(X[, 1]) - 100, max(X[, 1]) + 100)
ylim <- c(min(X[, 2]) - 100, max(X[, 2]) + 100)

# For each camera, get "time since first deployment began"
camera_locs <- camera_locs %>% 
  group_by(Date) %>% 
  mutate(mins_since_start = as.numeric(difftime(Time, min(Time), units = "mins")) + 20)

Dates <- unique(camera_locs$Date)

#### Data reformatting for model ####

# Make a matrix of observed counts
mtx_nrow <- max(table(camera_locs$Date))
ncam_vec <- numeric(3)
y_mtx <- matrix(ncol = 3, nrow = mtx_nrow)
for (i in 1:3) {
  ncam_vec[i] <- sum(camera_locs$Date == Dates[i])
  y_mtx[1:ncam_vec[i], i] <- camera_locs$count[camera_locs$Date == Dates[i]]
}

# Make a 3D array of camera coords  
X_array <- array(dim = c(mtx_nrow, 3, 2))
for (i in 1:3) {
  X_array[1:ncam_vec[i], i, 1] <- X$x[camera_locs$Date == Dates[i]]
  X_array[1:ncam_vec[i], i, 2] <- X$y[camera_locs$Date == Dates[i]]
}



data_template <- list(
  constants = list(M = M,
                   ncam = ncam_vec,
                   X = X_array,
                   xlim = xlim,
                   ylim = ylim),
  data = list(
    y = y_mtx
  ),
  inits = list(
    z = rbinom(M, 1, 0.5),
    psi = 0.5,
    p0 = 1,
    log_sigma = 3.5,
    s = array(0, dim = c(M, 3, 2))
  )
)

saveRDS(data_template, "simulations/data_template.RDS")
