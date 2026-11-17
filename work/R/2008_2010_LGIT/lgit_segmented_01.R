# breakpoints of influx
library(segmented)
year <- seq(10000,15000, by=50)
year2 <- year^2

influx <- read.csv("e:/Projects/igbpfire/yd/data/files/yd_influx_zts_b.csv")
attach(influx)
plot(ztinflux ~ est_age)

# basic trend
influx.lm1 <- lm(ztinflux ~ est_age, data=influx)
summary(influx.lm1)
plot(ztinflux ~ est_age)
lines(predict(influx.lm1, data.frame(est_age=year)) ~ year, col="red")

# segmented trend
influx.seg <- segmented(influx.lm1, seg.Z = ~ est_age, psi=list(est_age=c(11000,13000)), 
	control=seg.control(display=T, h=0.2) )
summary(influx.seg)
slope(influx.seg)

plot(ztinflux ~ est_age, pch=16, cex=0.5, xlim=c(15000,10000))
#plot(influx.seg, add=T, col="purple", lwd=3)
lines(influx.seg$fitted.values ~ est_age, col="purple", lwd=3)
lines(influx.seg, col="purple", lwd=2)

