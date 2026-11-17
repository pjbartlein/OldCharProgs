# correlations between PC scores and regional composite curves

pcapath <-"e:\\Projects\\GPWG\\analysis\\output\\v2n\\v2n_pcafiles\\"
loadname <- "z15-10_0100_50_ppca_04_load.csv" # "z22-0nrs_0500_50_ppca_04_load.csv"

loadfile <- paste(pcapath,loadname,sep="")

load <- read.csv(loadfile)

latfac <- cut(load$Lat, breaks=seq(-90,90,by=10))

plot(latfac, load$pc4, xlab="Latitude", xlim=c(1,18), ylab = "PC 1", ylim=c(-.25,.25), main=loadname)
abline(0,0)

plot(latfac, load$pc2, xlab="Latitude", xlim=c(1,18), ylab = "PC 2", ylim=c(-.25,.25), main=loadname)
abline(0,0)

plot(as.numeric(latfac), load$pc1, xlab="Latitude", xlim=c(1,18), ylab = "PC 1", ylim=c(-.25,.25), main=loadname, axes=F)
abline(0,0)
axis(2,seq(-.25,.25,by=.10), ylab="PC 1")
par(new=T)
plot(load$Lat, load$pc1, xlab=" ", ylab=" ", xlim=c(-90,90), axes=F, type="n")
axis(1,seq(-90,90,by=10))

plot(as.numeric(latfac), load$pc2, xlab="Latitude", xlim=c(1,18), ylab = "PC 2", ylim=c(-.25,.25), main=loadname, axes=F)
abline(0,0)
axis(2,seq(-.25,.25,by=.10), ylab="PC 2")
par(new=T)
plot(load$Lat, load$pc2, xlab=" ", ylab=" ", xlim=c(-90,90), axes=F, type="n")
axis(1,seq(-90,90,by=10))