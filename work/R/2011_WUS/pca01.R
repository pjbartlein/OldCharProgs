# pca of influx, peaks, tann and dai

workpath <-"e:\\Projects\\GPWG\\work\\wus\\R\\"

# interpolate peaks series

rawpeaksfile <- "peaks50_raw.csv"
outname <- "peaks50.csv"
outfile <- paste(workpath, outname, sep="")
                 
peaksin <- read.csv(paste(workpath, rawpeaksfile, sep=""))

targyear <- seq(from=500, to=2000, by=10)
peaks50 <- approx(peaksin$Year, peaksin$peak_density_050, targyear)
plot(peaksin$Year, peaksin$peak_density_050)
points(targyear, peaks50$y, pch=16, col="red")

peaksout <- cbind(targyear, peaks50$y)
colnames(peaksout) <- c("Year","Peaks50")
write.table(peaksout, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")

# pca

wusfile <- "combined.csv"
outname <- "pca.csv"
outfile <- paste(workpath, outname, sep="")

wusallin <- read.csv(paste(workpath, wusfile, sep=""))
summary(wusallin)

wusall <- wusallin[wusallin$Year <= 1800,]
wus <- wusall[,c(2,3,5,6)]

wus.pca <- princomp(wus, cor=T)
wus.pca
summary(wus.pca)
loadings(wus.pca)
screeplot(wus.pca)
plot(loadings(wus.pca))
biplot(wus.pca)