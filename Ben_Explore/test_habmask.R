library(terra)
library(tidyverse)
library(nimble)
library(nimbleSCR)


dMyHabMask <- nimbleFunction(run = function(x = double(0),
                                            s = double(1),
                                            xmax = double(0),
                                            xmin = double(0),
                                            ymax = double(0),
                                            ymin = double(0),
                                            resoln = double(0),
                                            habitatMask = double(2),
                                            log = logical(0)
                                            ){
  if (s[1] < xmin) return(-Inf)
  if (s[1] > xmax) return(-Inf)
  if (s[2] < ymin) return(-Inf)
  if (s[2] > ymax) return(-Inf)
  
  # Convert s to row/col
  row <- trunc((s[1] - xmin) / resoln) + 1
  col <- trunc((s[2] - ymin) / resoln) + 1
  
  test <- 1 - (habitatMask[row, col] == 0)
  
  if (log) {
    return(log(test)) 
  } else {
    return(test)
  }
  returnType(double(0))
})

rMyHabMask <- nimbleFunction(run = function(n = double(0),
                                            s = double(1),
                                            xmax = double(0),
                                            xmin = double(0),
                                            ymax = double(0),
                                            ymin = double(0),
                                            resoln = double(0),
                                            habitatMask = double(2)
                                            ){
  returnType(double(0))
  
  if (s[1] < xmin) return(0)
  if (s[1] > xmax) return(0)
  if (s[2] < ymin) return(0)
  if (s[2] > ymax) return(0)
  
  # Convert s to row/col
  row <- trunc((s[1] - xmin) / resoln) + 1
  col <- trunc((s[2] - ymin) / resoln) + 1
  
  if (habitatMask[row, col] == 0)  {
    return(0)
  } else {
    return(1)
  }
})
  
hb_raw <- rast("Data/Chicken_Rock_Map/ChickenRock_Classification.tif")

terra::values(hb_raw) <- ifelse(terra::values(hb_raw) == 4, 1, 0)
terra::values(hb_raw)[is.na(terra::values(hb_raw))] <- 0

hb_mask <- hb_raw %>% 
  aggregate(30, fun = "max")
terra::values(hb_mask)[is.nan(terra::values(hb_mask))] <- 0


hb_mtx <- t(as.matrix(hb_mask, wide = T))
hb_mtx <- hb_mtx[, ncol(hb_mtx):1]
hb_bbox <- ext(hb_mask)
hb_bbox <- as.numeric(c(hb_bbox[1], hb_bbox[2], hb_bbox[3], hb_bbox[4]))

hb_bbox[1:2] <- hb_bbox[1:2] - mean(hb_bbox[1:2])
hb_bbox[3:4] <- hb_bbox[3:4] - mean(hb_bbox[3:4])

testCode <- nimbleCode({
  for (i in 1:100) {
    s[i, 1] ~ dunif(xmin, xmax)
    s[i, 2] ~ dunif(ymin, ymax)
    ones[i] ~ dMyHabMask(
      s = s[i, 1:2], xmax = xmax, xmin = xmin, ymax = ymax, ymin = ymin,
      resoln = resoln, habitatMask = hb[1:nxcell, 1:nycell]
    )
  }
})


testMod <- nimbleModel(
  testCode,
  constants = list(
    xmin = hb_bbox[1], 
    xmax = hb_bbox[2], 
    ymin = hb_bbox[3], 
    ymax = hb_bbox[4],
    resoln = res(hb_mask)[1],
    hb = hb_mtx,
    nxcell = nrow(hb_mtx), 
    nycell = ncol(hb_mtx)
  ),
  data = list(
    ones = rep(1, 100)
  ),
  inits = list(s = cbind(
    rep(hb_bbox[1] + 30 * 21 + 1, 100),
    rep(hb_bbox[3] + 30 * 21 + 1, 100)
  )) # start in row 21 col 21
)
testMod$calculate()

cmod <- compileNimble(testMod)

conf <- configureMCMC(testMod, nodes = NA)

for (i in 1:100) {
  conf$addSampler(target = paste0("s[", i, ",1:2]"), type = "AF_slice")
}

mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc)

test_samples <- runMCMC(cmcmc, niter = 10000, nburnin = 0, thin = 1,
                        samplesAsCodaMCMC = TRUE)


plot(as.matrix(test_samples[, c("s[1, 1]", "s[1, 2]")]))

hb_bbox <- ext(hb_mask)
hb_bbox <- as.numeric(c(hb_bbox[1], hb_bbox[2], hb_bbox[3], hb_bbox[4]))

pts <-as.matrix(test_samples[, c("s[1, 1]", "s[1, 2]")]) %>% 
  as.data.frame() %>% 
  rename(x = `s[1, 1]`, y = `s[1, 2]`) %>% 
  mutate(x = x + mean(hb_bbox[1:2]), y = y + mean(hb_bbox[3:4])) %>% 
  vect(geom = c("x", "y"), crs = crs(hb_mask))

plot(hb_mask)
points(pts)

