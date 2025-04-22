process_ROV_buffers <- function(buffer_size, surface_to_buffer) {
  
  rov_locs <- read_csv("Data/SnapperMvmtAbundanceStudy/CountData/ROV/ROV_latlong_real_2023_processed.csv")
  rov_counts <- read_xlsx("Data/SnapperMvmtAbundanceStudy/CountData/ROV/RS_2023_2024_ROV_counts_NC_FL_VPS_studies.xlsx", skip = 1) %>% 
    filter(Region == 1)
  
  rov_dat <- bind_rows(
    rov_counts %>% 
      dplyr::select(
        Site_ID = `Site Number`,
        Date = `Date`,
        Latitude, 
        Longitude,
        depth = `ROV Depth (ft)`,
        FOV,
        # length = `Transect Length1`,
        area = `Transect 1 Area`,
        hab_quality = `Habitat Quality`,
        reef_pct = `Estimated % Transect over reef`,
        count = `Transect 1 RS count`
      ) %>% 
      mutate(transect = 1),
    rov_counts %>% 
      dplyr::select(
        Site_ID = `Site Number`,
        Date = `Date`,
        Latitude, 
        Longitude,
        depth = `ROV Depth (ft)`,
        FOV,
        # length = `Transect Length1`,
        area = `Transect 2 Area`,
        hab_quality = `Habitat Quality`,
        reef_pct = `Estimated % Transect over reef`,
        count = `Transect 2 RS count`
      ) %>% 
      mutate(transect = 2)
  ) %>% 
    mutate(density = count/area, 
           structured_area = area * reef_pct)
  
  
  
  
  df_split <- split(rov_locs, 
                    paste0(rov_locs$Site_ID, rov_locs$Date, rov_locs$transect))
  line_geoms <- lapply(df_split, 
                       function(group) {
                         as.matrix(group[, c("longitude", "latitude")])
                       })
  
  # Create the line SpatVector
  lines_vect <- vect(line_geoms, type = "lines", crs = "+proj=longlat") %>% 
    project(crs(surface_to_buffer))
  
  line_attributes <- do.call(rbind, lapply(df_split, function(group) {
    group[1, c("Site_ID", "Date", "transect")]  # retain one row per group
  }))
  
  line_attributes$buffer_area <- NA
  line_attributes$surface_mean <- NA
  line_attributes$surface_sum <- NA
  
  intersections_df <- data.frame()
  
  for (i in 1:nrow(line_attributes)) {
    this_line <- lines_vect[i, ]
    
    buffer <- buffer(this_line, buffer_size)
    
    this_value <- terra::extract(surface_to_buffer, buffer, weights = TRUE, cells = TRUE)
    
    line_attributes$buffer_area[i]  <- expanse(buffer)
    line_attributes$surface_mean[i] <- sum(this_value$lyr.1 * (this_value$weight / sum(this_value$weight)))
    line_attributes$surface_sum[i]  <- sum(this_value$lyr.1 * this_value$weight)
    
    this_value <- bind_cols(this_value, line_attributes[i, c("Site_ID", "Date", "transect")])
    intersections_df <- bind_rows(intersections_df, this_value)
  }
  
  intersections_df <- intersections_df %>% 
    mutate(
      y_ind = 1 + dim(surface_to_buffer)[1] - ceiling(cell / dim(surface_to_buffer)[2]),
      x_ind = (cell %% dim(surface_to_buffer)[2]),
      ROV_ID = paste(Site_ID, Date, transect, sep = "_")
    ) %>% 
    arrange(ROV_ID)
  
  
  rov_dat <- rov_dat %>% 
    left_join(line_attributes, by = c("Site_ID", "Date", "transect")) %>% 
    mutate(ROV_ID = paste(Site_ID, Date, transect, sep = "_")) %>% 
    arrange(ROV_ID)
  
  rbs <- rbe <- numeric(nrow(rov_dat))
  for (i in 1:length(rbs)) {
    rbs[i] <- min(which(intersections_df$ROV_ID == rov_dat$ROV_ID[i]))
    rbe[i] <- max(which(intersections_df$ROV_ID == rov_dat$ROV_ID[i]))
  }
  
  return(list(
    rbe = rbe,
    rbs = rbs,
    rov_dat = rov_dat,
    intersections_df = intersections_df
  ))
}


