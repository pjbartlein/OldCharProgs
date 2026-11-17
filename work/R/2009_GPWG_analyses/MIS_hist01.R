# histograms

# read data e:\Projects\GPWG\data\files\v2m\v2m-lr_all_influx2.csv
v2m.lr <- read.csv(file.choose())
attach(v2m.lr)

# set missing data values to NA
tall[tall == -9999] <- NA
ztall[ztall == -9999] <- NA

oldpar <- par(mfcol=c(3,2))
hist(ztall[est_age > -60 & est_age <= 11700], breaks=seq(-10,+10,by=0.2), xlim=c(-10,10), main="Holocene")
hist(ztall[est_age > 11700 & est_age <= 73500], breaks=seq(-10,+10,by=0.2), xlim=c(-10,10), main="Glacial")
hist(ztall[est_age > 115000 & est_age <= 128000], breaks=seq(-10,+10,by=0.2), xlim=c(-10,10), main="Eemian")

hist(ztall[est_age > 14700 & est_age <= 27800], breaks=seq(-10,+10,by=0.2), xlim=c(-10,10), main="MIS 2")
hist(ztall[est_age > 27800 & est_age <= 59400], breaks=seq(-10,+10,by=0.2), xlim=c(-10,10), main="MIS 3")
hist(ztall[est_age > 59400 & est_age <= 73500], breaks=seq(-10,+10,by=0.2), xlim=c(-10,10), main="MIS 4")

oldpar <- par(mfcol=c(3,2))
hist(tall[est_age > -60 & est_age <= 11700], breaks=seq(0,1,by=0.025), xlim=c(0,1), main="Holocene")
hist(tall[est_age > 11700 & est_age <= 73500], breaks=seq(0,1,by=0.025), xlim=c(0,1), main="Glacial")
hist(tall[est_age > 115000 & est_age <= 128000], breaks=seq(0,1,by=0.025), xlim=c(0,1), main="Eemian")

hist(tall[est_age > 14700 & est_age <= 27800], breaks=seq(0,1,by=0.025), xlim=c(0,1), main="MIS 2")
hist(tall[est_age > 27800 & est_age <= 59400], breaks=seq(0,1,by=0.025), xlim=c(0,1), main="MIS 3")
hist(tall[est_age > 59400 & est_age <= 73500], breaks=seq(0,1,by=0.025), xlim=c(0,1), main="MIS 4")
