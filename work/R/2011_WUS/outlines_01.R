# draw "squarish" coastlines given a gridded land/ocean/ice mask

library(ncdf)

# check and reset working directory if necessary
getwd()
#workdir <- "e:\\DataVis\\raster\\"
#setwd(workdir)

ncfname <- "21ka_points_ice5g-dqs_fixed_mask.nc"
dname <- "landmask"

landoutlinefile <- "21ka_land_hd.csv" 
iceoutlinefile <- "21ka_ice_hd.csv" 

# open the netCDF file
ncin <- open.ncdf(ncfname) 

# print some basic informaition
print(ncin)

# more verbose description
# str(ncin)

# get longitudes and latitudes (grid-cell centers)
lon <- get.var.ncdf(ncin,"lon")
nlon <- dim(lon)
head(lon)
lat <- get.var.ncdf(ncin,"lat")
nlat <- dim(lat)
head(lat)
print(c(nlon,nlat))

# get the mask and attributes
landmask <- get.var.ncdf(ncin,dname)
dlname <- att.get.ncdf(ncin,dname,"long_name")
code <- att.get.ncdf(ncin,dname,"code")
dim(landmask)

# done with the netCDF file, so close it
close.ncdf(ncin)

# get lats and lons of cell edges
lon.e <- dim(nlon); lon.w <- dim(nlon)
for (j in 1:nlon-1) {
  lon.e[j] <- (lon[j]+lon[j+1])/2.0
  lon.w[j+1] <- (lon.e[j])
}
lon.w[1] <- -180.0
lon.e[nlon] <- 180.0
head(cbind(lon.w,lon,lon.e))
tail(cbind(lon.w,lon,lon.e))

lat.n <- dim(nlat); lat.s <- dim(nlat)
for (k in 1:nlat-1) {
  lat.n[k] <- (lat[k]+lat[k+1])/2.0
  lat.s[k+1] <- (lat.n[k])
}
lat.s[1] <- -90.0
lat.n[nlat] <- 90.0
head(cbind(lat.s,lat,lat.n))
tail(cbind(lat.s,lat,lat.n))

# get land and ice points
land <- landmask
land[land==2] <- 0
ice <- landmask
ice[ice==1] <- 0

# loop over lats and lons and get land outline segments
# latitudes
xy <- matrix(nrow=nlon*nlat*2,ncol=4)
i <- 0
for (k in 1:nlat) {
  for (j in 1:(nlon-1)) {  
    if(land[j,k] != land[j+1,k]) {
      i <- i+1
      xy[i,1] <- lon.e[j]; xy[i,2] <- lat.s[k]; xy[i,3] <- lon.e[j]; xy[i,4] <- lat.n[k]
    }
  }
}

# longitudes
for (j in 1:nlon) {
  for (k in 1:(nlat-1)) {  
    if(land[j,k] != land[j,k+1]) {
      i <- i+1
      xy[i,1] <- lon.w[j]; xy[i,2] <- lat.n[k]; xy[i,3] <- lon.e[j]; xy[i,4] <- lat.n[k]
    }
  }
}
land.segs <- data.frame(na.omit(xy))
names(land.segs) <- c("x0", "y0", "x1", "y1")
head(land.segs)
tail(land.segs)

# now ice

# loop over lats and lons and get land outline segments
# latitudes
xy <- matrix(nrow=nlon*nlat*2,ncol=4)
i <- 0
for (k in 1:nlat) {
  for (j in 1:(nlon-1)) {  
    if(ice[j,k] != ice[j+1,k]) {
      i <- i+1
      xy[i,1] <- lon.e[j]; xy[i,2] <- lat.s[k]; xy[i,3] <- lon.e[j]; xy[i,4] <- lat.n[k]
    }
  }
}

# longitudes
for (j in 1:nlon) {
  for (k in 1:(nlat-1)) {  
    if(ice[j,k] != ice[j,k+1]) {
      i <- i+1
      xy[i,1] <- lon.w[j]; xy[i,2] <- lat.n[k]; xy[i,3] <- lon.e[j]; xy[i,4] <- lat.n[k]
    }
  }
}
ice.segs <- data.frame(na.omit(xy))
names(ice.segs) <- c("x0", "y0", "x1", "y1")
head(ice.segs)
tail(ice.segs)

plot(NULL,NULL,xlim=c(-180,180),ylim=c(-90,90))
segments(land.segs$x0,land.segs$y0,land.segs$x1,land.segs$y1, lwd=2)
segments(ice.segs$x0,ice.segs$y0,ice.segs$x1,ice.segs$y1, col="magenta", lwd=2)

# save segments for reuse
write.table(land.segs, landoutlinefile, row.names=FALSE, sep=",")
write.table(ice.segs, iceoutlinefile, row.names=FALSE, sep=",")

