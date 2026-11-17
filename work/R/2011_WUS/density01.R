pfile <- read.csv("e:\\Projects\\GPWG\\data\\peaks\\WUS\\files\\wus_peaks_ts.csv")
attach(pfile)
summary(pfile)

bw=25.0
begyr=-50.0; endyr=3000.0
hist(age[peaks==1], freq=F,  breaks=seq(begyr-(2*bw),endyr+2*bw,2*bw), 
    xlim=c(endyr+bw,begyr-bw), main=paste(" bw =", as.character(bw), "\n"))


lines(density(age[[peaks==1], bw=bw), col="red")




cbind(seq(1,length(age[finalpks==1])),site[finalpks==1],age[finalpks==1])