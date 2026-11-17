# piecewise regression fits
library(segmented)

workpath <-"e:\\Projects\\GPWG\\analysis\\output\\v2n\\v2n_curves2\\"
curvesname <- "v2n_lgit_0100_b.csv"

infile <- paste(workpath,curvesname,sep="")
outfile <- paste(workpath,"v2n_lgit_0100_b_pwfits.csv",sep="")

lgit <- read.csv(infile)
lgit <- lgit[lgit$targ_age >= 10000.0 & lgit$targ_age <= 15000.0,]
year <- seq(10000,15000, by=20)
attach(lgit)

# basic trend
lm_na_1 <- lm(na ~ targ_age, data=lgit)
summary(lm_na_1)
plot(na ~ targ_age, xlim=c(15000,10000), pch=16)
lines(predict(lm_na_1, data.frame(targ_age=year)) ~ year, lwd=3, col="red")

# segmented trend
seg_na_1 <- segmented(lm_na_1, seg.Z = ~ targ_age, psi=list(targ_age=c(11000,13000)), 
	control=seg.control(display=T, h=0.1) )
summary(seg_na_1)
slope(seg_na_1)

plot(ztinflux ~ targ_age, pch=16, cex=0.5, xlim=c(15000,10000))
#plot(influx.seg, add=T, col="purple", lwd=3)
lines(influx.seg$fitted.values ~ targ_age, col="purple", lwd=3)
lines(influx.seg, col="purple", lwd=2)

