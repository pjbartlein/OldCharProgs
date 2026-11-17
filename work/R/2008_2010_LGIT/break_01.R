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

# segmented polynomial trend 2
influx.lm6 <- lm(ztinflux ~ est_age + I(est_age^2), data=influx)
summary(influx.lm6)
plot(ztinflux ~ est_age, xlim=c(15000,10000)
lines(predict(influx.lm6, data.frame(est_age=year)) ~ year, col="red", lwd=2)

influx.seg2 <- segmented(influx.lm6, seg.Z = ~ est_age, psi=list(est_age=c(11600,14000)), 
	control=seg.control(display=T, h=0.2) )
summary(influx.seg2)
slope(influx.seg2)

plot(ztinflux ~ est_age, pch=16, cex=0.5, xlim=c(15000,10000))
#plot(influx.seg2, add=T, col="purple", lwd=3)
lines(influx.seg2$fitted.values ~ est_age, col="purple", lwd=2)
lines(influx.seg2, col="purple", lwd=2)

# segmented polynomial trend 3
influx.lm7 <- lm(ztinflux ~ est_age + I(est_age^2) + I(est_age^3), data=influx)
summary(influx.lm7)
#plot(ztinflux ~ est_age)
lines(predict(influx.lm7, data.frame(est_age=year)) ~ year, col="black", lwd=2)

influx.seg3 <- segmented(influx.lm7, seg.Z = ~ est_age, psi=list(est_age=c(11700,12900)), 
	control=seg.control(display=T, h=0.2) )
summary(influx.seg3)
slope(influx.seg3)

plot(ztinflux ~ est_age, pch=16, cex=0.5, xlim=c(15000,10000))
lines(influx.seg3$fitted.values ~ est_age, col="blue", lwd=2)
lines(influx.seg3, col="blue", lwd=2)