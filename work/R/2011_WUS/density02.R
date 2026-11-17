pfile <- read.csv("e:\\Projects\\GPWG\\work\\wus\\data\\files\\wus_peaks_ts.csv")

attach(pfile)
summary(pfile)

bw=50.0
begyr=-75.0; endyr=3200.0
hist(age[peaks==1], freq=F,  breaks=seq(begyr-(2*bw),endyr+2*bw,2*bw), 
    xlim=c(endyr+bw,begyr-bw), main=paste(" bw =", as.character(bw), "\n"))

begyr=-50.0; endyr=3000.0
peakden <- density(age[peaks==1], bw=bw, from=begyr, to=endyr)
lines(peakden, xlim=c(endyr-bw,begyr+bw), col="blue", lwd=2)

peakdenout <- cbind(peakden$x,peakden$y)
colnames(peakdenout) <- c("Year","Density")

outfile <- "e:\\Projects\\GPWG\\work\\wus\\data\\files\\wus_peaks_050.csv"
write.table(peakdenout, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")
