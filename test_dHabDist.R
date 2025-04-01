library(nimbleEcology)
library(tidyverse)
library(readxl)
library(terra)

dHabDistr <- nimbleFunction(run = function(x = double(1),
                                           xmax = double(0),
                                           xmin = double(0),
                                           ymax = double(0),
                                           ymin = double(0),
                                           resoln = double(0),
                                           habitatMask = double(2),
                                           log = logical(0)) {
  if (x[1] < xmin) return(-Inf)
  if (x[1] > xmax) return(-Inf)
  if (x[2] < ymin) return(-Inf)
  if (x[2] > ymax) return(-Inf)
  
  # Convert s to row/col
  row <- trunc((x[1] - xmin) / resoln) + 1
  col <- trunc((x[2] - ymin) / resoln) + 1
  
  if (log) return(log(habitatMask[row, col]))
  else return(habitatMask[row, col])
  returnType(double(0))
}, buildDerivs = FALSE)

rHabDistr <- nimbleFunction(run = function(n = double(0),
                                           xmax = double(0),
                                           xmin = double(0),
                                           ymax = double(0),
                                           ymin = double(0),
                                           resoln = double(0),
                                           habitatMask = double(2)
){
  returnType(double(1))
  
  # Randomly select a row based on their relative weights
  rowProbs <- numeric(dim(habitatMask)[1])
  for (i in 1:(dim(habitatMask)[1])) {
    rowProbs[i] <- sum(habitatMask[i, ])
  }
  row <- rcat(1, prob = rowProbs)
  
  # Randomly select a column from that row
  col <- rcat(1, prob = habitatMask[row, ])
  
  # # Generate a random coord. uniformly within that cell
  x <- c(
    xmin + resoln * (row - 1) + runif(1, 0, resoln),
    ymin + resoln * (col - 1) + runif(1, 0, resoln)
  )
  return(x)
  
}, buildDerivs = FALSE)
# ctest <- compileNimble(rHabDistr)

initialize_s <- function(n, 
                         xmax,
                         xmin,
                         ymax,
                         ymin,
                         resoln,
                         habitatMask) {
  s <- matrix(NA, nrow = n, ncol = 2)
  for (i in 1:n) {
    s[i, ] <- rHabDistr(1, xmax = xmax, xmin = xmin,
                        ymax = ymax, ymin = ymin, resoln = resoln,
                        habitatMask = habitatMask)
  }
  s
}

testcode <- nimbleCode({
  for (i in 1:10) {
    x[i, 1:2] ~ dHabDistr(
      xmax = 1,
      xmin = 0,
      ymax = 1,
      ymin = 0,
      resoln = 0.1,
      habitatMask = hm[1:10, 1:10]
    )
  }
})

test_mtx <- matrix(runif(100, 0, 1), nrow = 10)
test_mtx <- test_mtx/sum(test_mtx)

testmod <- nimbleModel(testcode, 
                       constants = list(hm = test_mtx),
                       inits = list(
                         x = matrix(runif(20, 0, 1), nrow = 10)
                       ))

cmod <- compileNimble(testmod)

mcmcConf <- configureMCMC(testmod)
mcmc <- buildMCMC(mcmcConf)
cmcmc <- compileNimble(mcmc)

samples <- runMCMC(cmcmc, niter = 100000, nburnin = 0, thin = 1,
                   samplesAsCodaMCMC = TRUE)

plot(as.numeric(samples[, "x[1, 1]"]), as.numeric(samples[, "x[1, 2]"]), cex = 0.2)
image(test_mtx)


#### Now test with the real distribution ####

hb_raw <- rast("Data/Chicken_Rock_Map/ChickenRock_Classification.tif")

# Get all the VPS fixes of living fish
VPS_folder <- "Data/SnapperMvmtAbundanceStudy/VPS_Data/VPS-ChickenRock-01-Results-20240202/results/animal/"
VPS_files <- list.files(VPS_folder)
fish_files <- VPS_files[grepl("^12", VPS_files)]
fish_fate <- read_xlsx("Data/SnapperMvmtAbundanceStudy/VPS_Data/Fate_assignments_from_discard_paper_vps 2023.xlsx")
fish_to_drop <- fish_fate$`Tag number`[fish_fate$`subsequent assignments based on velocity and depth data` == "Discard mortality"]
dead_fish_files <- paste0(fish_to_drop, ".csv")
fish_files <- fish_files[!fish_files %in% dead_fish_files]
fish_positions <- lapply(file.path(VPS_folder, fish_files),
                         read_csv) %>% 
  bind_rows()
fish_pts <- vect(fish_positions, geom = c("Longitude", "Latitude"),
                 crs = "+proj=longlat") %>% 
  project(crs(hb_raw))

e <- ext(fish_pts)
xmin(e) <- floor(xmin(e) / 50) * 50
xmax(e) <- ceiling(xmax(e) / 50) * 50
ymin(e) <- floor(ymin(e) / 50) * 50
ymax(e) <- ceiling(ymax(e) / 50) * 50

template_grid <- rast(e, res = 50)
terra::values(template_grid) <- 1:ncell(template_grid)


cell_counts <- count(extract(template_grid, fish_pts), lyr.1)
cell_counts <- left_join(data.frame(lyr.1 = 1:ncell(template_grid)), cell_counts)
cell_counts$n[is.na(cell_counts$n)] <- 0
cell_counts$prob <- (cell_counts$n+1) / sum(cell_counts$n)

vps_intensity_ras <- template_grid
terra::values(vps_intensity_ras) <- cell_counts$prob
plot(vps_intensity_ras)

vps_mtx <- t(as.matrix(vps_intensity_ras, wide = T))
vps_mtx <- vps_mtx[, ncol(vps_mtx):1]

grid_bbox <- ext(vps_intensity_ras)
grid_bbox <- as.numeric(c(grid_bbox[1], grid_bbox[2], grid_bbox[3], grid_bbox[4]))

x_offset <- mean(grid_bbox[1:2])
y_offset <- mean(grid_bbox[3:4])

grid_bbox[1:2] <- grid_bbox[1:2] - x_offset
grid_bbox[3:4] <- grid_bbox[3:4] - y_offset



testCode <- nimbleCode({
  for (i in 1:10) {
    x[i, 1:2] ~ dHabDistr(
      xmax = xmax,
      xmin = xmin,
      ymax = ymax,
      ymin = ymin,
      resoln = resoln,
      habitatMask = hm[1:hm_nrow, 1:hm_ncol]
    )
  }
})

mod <- nimbleModel(
  testCode,
  constants = list(
    xmin = grid_bbox[1], 
    xmax = grid_bbox[2], 
    ymin = grid_bbox[3], 
    ymax = grid_bbox[4],
    resoln = 50, 
    hm = vps_mtx,
    hm_nrow = nrow(vps_mtx),
    hm_ncol = ncol(vps_mtx)
  ),
  inits = list(
    x = initialize_s(n = 10,
                     xmin = grid_bbox[1],
                     xmax = grid_bbox[2],
                     ymin = grid_bbox[3],
                     ymax = grid_bbox[4],
                     resoln = 50, habitatMask = vps_mtx)
  )
)


cmod <- compileNimble(mod)

mcmcConf <- configureMCMC(mod)
mcmc <- buildMCMC(mcmcConf)
cmcmc <- compileNimble(mcmc)

samples <- runMCMC(cmcmc, niter = 100000, nburnin = 0, thin = 1,
                   samplesAsCodaMCMC = TRUE)

plot(as.numeric(samples[, "x[1, 1]"]), as.numeric(samples[, "x[1, 2]"]), cex = 0.2)
image(vps_mtx)
