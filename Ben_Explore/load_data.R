rm(list=ls())
library(class)
library(raster)
library(fields)
library(sf)
library(RANN)

# Plotting functions

quilt_plot_atl  <- function(s,y,main=NULL,xlim=c(-83.97,-74.74),ylim=c(24.33,36.54)){
  library(maps)
  quilt.plot(s[,1],s[,2],y,nx=100,ny=100,asp=1,
             xlim=xlim,ylim=ylim,col=viridis(50),
             main=main,xlab="",ylab="")
  map("state",add=TRUE)
}

plot_atl  <- function(s,y,main=NULL,xlim=c(-83.97,-74.74),ylim=c(24.33,36.54)){
  plot(s,col=y,asp=1,pch=19,xlim=xlim,ylim=ylim,main=main,xlab="",ylab="")
  map("state",add=TRUE)
}


# Load the study area and raster

study_area   <- st_read("Data/StudyArea/StudyArea.shp")
rast         <- raster("Data/Hardbottom/ProcessedRasters/r_TNC.tif")
values(rast) <- ifelse(is.na(values(rast)),0,values(rast))
rast         <- mask(rast,study_area)

# Load the TNC hard bottom data and distance to HB

keep       <- which(!is.na(values(rast)))
grid_hb    <- values(rast)[keep]
grid_s     <- coordinates(rast)[keep,]
grid_area  <- area(rast)[keep]
grid_d2hb  <- nn2(grid_s[grid_hb>0,], query = grid_s,k=1)$nn.dists
rm(keep,rast)

# Load the depth data

dp     <- raster("Data/depth/depth_wgs.tif")
depth  <- extract(dp, grid_s, method="bilinear") 
isna   <- which(is.na(depth))
notna  <- which(!is.na(depth)) 
knn    <- knn(train = coordinates(dp)[notna,],
              test  = coordinates(dp)[isna,],
              cl    = depth[notna],k=1)
depth[isna] <- as.numeric(as.character(knn)) 
grid_depth  <- depth
rm(dp,knn,isna,notna,depth)


# ROV count data

rov_dat    <- read.csv("Data/Counts/cleanROVdata2021_2022.csv")
rov_ID     <- rov_dat$SiteID
rov_s      <- rov_dat[,c("Longitude", "Latitude")]
rov_year   <- rov_dat$Year
rov_cover1 <- rov_dat$T1cover
rov_cover2 <- rov_dat$T2cover
rov_area   <- rowSums(rov_dat[,c("T1AreaSurv","T2AreaSurv")])/(1000^2)
rov_count  <- rowSums(rov_dat[,c("T1Count","T2Count")])
rov_cell   <- NULL 
rov_rand   <- rov_dat$RandomSelection=="Y"
rov_cell   <- nn2(grid_s, query = rov_s, k=1)$nn.idx
rov_hb     <- grid_hb[rov_cell]
rov_d2hb   <- grid_d2hb[rov_cell,]
rov_depth  <- grid_depth[rov_cell]

# Camera count data

cam_dat     <- read.csv("Data/Counts/cleanCameraData2011_2022.csv")
cam_dat     <- cam_dat[cam_dat$Year>2020,]
cam_ID      <- cam_dat$SiteID
cam_s       <- cam_dat[,c("Start_Longitude", "Start_Latitude")]
cam_subst   <- cam_dat$Substrate
cam_year    <- cam_dat$Year
cam_count   <- as.matrix(cam_dat[,20:60])
cam_mean    <- apply(cam_dat[,19 + 1:41],1,mean)
cam_median  <- apply(cam_dat[,19 + 1:41],1,median)
cam_max     <- apply(cam_dat[,19 + 1:41],1,max)
cam_cell    <- nn2(grid_s, query = cam_s, k=1)$nn.idx
cam_cur_dir <- cam_dat$Current_Direction
cam_turb    <- cam_dat$Turbidity+1
cam_area    <- rep(0,nrow(cam_s))
cam_area[which(cam_dat$Current_Direction==0)] <- 0.048
cam_area[which(cam_dat$Current_Direction==1)] <- 0.057
cam_area[which(cam_dat$Current_Direction==2)] <- 0.007
cam_hb      <- grid_hb[cam_cell]
cam_d2hb    <- grid_d2hb[cam_cell,]
cam_depth   <- grid_depth[cam_cell]

# Compute CAR basis functions

thin <- 10000
v    <- grid_s[1:nrow(grid_s)%%thin==0,]
nv   <- nrow(v)
nn   <- nn2(v,query=v,k=5)$nn.idx[,-1]
A    <- matrix(0,nv,nv)
for(i in 1:nv){
  A[i,nn[i,]] <- 1
  A[nn[i,],i] <- 1
}
CAR      <- as.spam(diag(rowSums(A))-A)
CAR      <- eigen(CAR)
CAR      <- sweep(CAR$vec[,nv-1:100-1],2,sqrt(CAR$val[nv-1:100-1]),"/")
CAR_cell <- nn2(v,query=grid_s,k=1)$nn.idx
colnames(CAR) <- paste0("CARbasis",1:100)
rm(v,nv,nn,A,i,thin)

# Export combined data

filename <- "Data/combined_data3.RData"
save.image(file=filename)

print("Crude estimate from ROV")
print(sum(grid_area)*sum(rov_count)/sum(rov_area))

print("Crude estimate from ROV - random sites")
print(sum(grid_area)*sum(rov_count[rov_rand])/sum(rov_area[rov_rand]))

print("Crude estimate from ROV - nonrandom sites")
print(sum(grid_area)*sum(rov_count[!rov_rand])/sum(rov_area[!rov_rand]))

print("Crude estimates (mean/median/max) from Camera")
print(sum(grid_area)*sum(cam_mean)/sum(cam_area))
print(sum(grid_area)*sum(cam_median)/sum(cam_area))
print(sum(grid_area)*sum(cam_max)/sum(cam_area))
