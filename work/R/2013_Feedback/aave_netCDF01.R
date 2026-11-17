## netCDF

library(ncdf)

# calculate area-weighted averages of time series from an nlon x nlat x nt
# netCDF file

# check and reset working directory if necessary
getwd()
#workdir <- "e:\\DataVis\\raster\\"
#setwd(workdir)

ncfname <- "Mann2009_tmp.nc"
dname <- "tmp_anm"

ptsfile <- "china_pts.csv" # lons and lats of points in the region
aavefile <- "aaves_china_pts_tmp_R.csv"
aavefile2 <- "aaves_china_pts_tmp_R_2.csv"

# open the netCDF file
ncin <- open.ncdf(ncfname) 

# print some basic informaition
print(ncin)

# more verbose description
# str(ncin)

# get longitudes and latitudes
lon <- get.var.ncdf(ncin,"lon")
nlon <- dim(lon)
head(lon)
lat <- get.var.ncdf(ncin,"lat")
nlat <- dim(lat)
head(lat)
print(c(nlon,nlat))

# get time (cft:  ClimateForecast convention for time)
cft <- get.var.ncdf(ncin,"time")
tunits <- att.get.ncdf(ncin,"time","units")
head(t,20); tail(cft,20)
nt <- dim(cft)
nt
tunits$value

# get areas of grid cells
area <- get.var.ncdf(ncin,"area")
head(area)
length(area)

# now get the data and attributes
tmp.array <- get.var.ncdf(ncin,dname)
dlname <- att.get.ncdf(ncin,dname,"long_name")
dunits <- att.get.ncdf(ncin,dname,"units")
fillvalue <- att.get.ncdf(ncin,dname,"_FillValue")
dim(tmp.array)

# get global attributes, if desired
title <- att.get.ncdf(ncin,0,"title")
institution <- att.get.ncdf(ncin,0,"institution")
datasource <- att.get.ncdf(ncin,0,"source")
citation <- att.get.ncdf(ncin,0,"citation")
history <- att.get.ncdf(ncin,0,"history")
Conventions <- att.get.ncdf(ncin,0,"Conventions")

# done with the netCDF file, so close it
close.ncdf(ncin)

# show "real" times
tustr <- strsplit(tunits$value, " ")
tdstr <- strsplit(unlist(tustr)[3], "-")
tyear=as.integer(unlist(tdstr)[1])
year <- cft+tyear
head(year); tail(year)

# replace netCDF fillvalues with R NAs
tmp.array[tmp.array==fillvalue$value] <- NA

# get target points that define the region
pts <- read.csv(ptsfile)
pts
npts <- dim(pts)[1]
npts


# one way to get area-weighted averages: leave data in the 3-d array
# and get area averages of individual time slices in a loop

# subscripts of points in the region
j <- match(pts$tlon, lon)
k <- match(pts$tlat, lat)
cbind(j,k)

# define array time-slice area-average function
aave.array <- function(tmp.in,j,k,n,area.in) {
  tmp.region <- tmp.in[cbind(j,k,n)]
  area.region <- area.in[cbind(j,k,n)[,2]]
  aave.array <- weighted.mean(tmp.region, area.region, na.rm=TRUE) 
}

# demonstration of aave.array() for first time step
#n <-1 
#tmp.region <- tmp.array[cbind(j,k,n)]
#area.region <- area[cbind(j,k,n)[,2]]
#aave.demo <- weighted.mean(tmp.region, area.region, na.rm=TRUE) 
#aave.demo

# get the time series of area-weighted averages
aave <- dim(nt)
for (n in 1:nt) {
  aave[n] <- aave.array(tmp.array,j,k,n,area)
}

plot(aave ~ year, type="l")

# write out the area averages
write.table(na.omit(cbind(year,aave)),aavefile, row.names=FALSE, sep=",")



# another approach -- make a nt rows x npts cols data frame

# convert the nlon by nlat by nt array into a nlon*lat (rows) by nt (cols) matrix
tmp.vec.long <- as.vector(tmp.array) # flatten the array
length(tmp.vec.long)
tmp.mat <- matrix(tmp.vec.long, nrow=nlon*nlat, ncol=nt)
dim(tmp.mat)

# lons and lats (and j's and k's of each row)
lonlat <- expand.grid(lon,lat)
dim(lonlat)
head(lonlat)
jk <- expand.grid(seq(1:nlon),seq(1:nlat))
dim(jk)
head(cbind(jk,lonlat))
tail(cbind(jk,lonlat))

# create a row-selection variable
jk2 <- (jk[,1]*100)+jk[,2]
head(jk2)
tail(jk2)

# rows to select (from j's and k's of target points)
jksel <- (j*100)+k
jksel

# test selection with lon and lats (will also be used or names)
rowselect <- match(jksel, jk2)
lonlatselect <- lonlat[rowselect,]
lonlatselect

# select the appropriate rows from tmp.mat
tmp.mat.select <- tmp.mat[rowselect,]
dim(tmp.mat.select)
tmp.mat.select[,1:10]

# transpose to make an nt (rows) by nlon*nlat (cols) matrix
tmp.mat.select <- t(tmp.mat.select)
dim(tmp.mat.select)
head(tmp.mat.select)

# grid-cell areas of the target points
areapts <- area[k]
areapts


# define data-frame row area-average function
# areapts is a global variable in this function
aave.row <- function(tmp.in) {
  aave.row <- weighted.mean(tmp.in, areapts, na.rm=TRUE) }

# demonstration of aave.row()
#n <- 1
#aave2.demo <-weighted.mean(tmp.mat.select[n,][1:npts],areapts, na.rm=TRUE)
#aave2.demo

# area averages of rows in data frame
aave2 <- dim(nt)
aave2 <- apply(tmp.mat.select,1,aave.row)
head(aave2)
plot(year, aave2, type="l", col="red")


# make a data frame
tmp.datfrm <- data.frame(cbind(year,tmp.mat.select,aave2))
head(tmp.datfrm)

# names
lonlatselect <- t(lonlatselect)
dim(lonlatselect)
lonlatselect
lonlatnames <- paste(as.character(lonlatselect[1,]),"E_",as.character(lonlatselect[2,]),"N",sep="")
lonlatnames

colnames(tmp.datfrm) <- c("year",lonlatnames,"areaave")
head(tmp.datfrm)
summary(tmp.datfrm)

# write out the data frame including area averages
write.table(tmp.datfrm, aavefile2, row.names=FALSE, sep=",")

