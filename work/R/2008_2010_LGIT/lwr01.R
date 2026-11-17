# lwr of YD influx data

workpath <- "e:\\Projects\\igbpfire\\yd\\data2\\"
infile <- "yd_influx_zts.csv"

zt <- read.csv(paste(workpath, infile, sep=""))
summary(zt)
attach(zt)

plot(ztinflux ~ est_age, pch=16, col="gray", cex=0.5, xlim=c(16000, 9000), ylim=c(-5,5))

new.age <- data.frame(est_age=seq(10000, 15000, by=20))

loess.model <- loess(ztinflux ~ est_age, span=0.05)
loess.model
zthat <- predict(loess.model, newdata=new.age)
plot(new.age$est_age, zthat, type="l", col="black", ylim=c(-2,3), xlim=c(15000, 10000),
     xlab="Age (cal yr BP)", ylab="Z-Scores of Transformed Charcoal Influx")
points(new.age$est_age, zthat, pch=16, cex=0.5, col="red")
