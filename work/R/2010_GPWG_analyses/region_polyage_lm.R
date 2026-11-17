# regressions of charcoal curves on PCs
workpath <-"e:\\Projects\\GPWG\\analysis\\output\\v2n\\v2n_curves2\\"
curvesname <- "v2n_all_nsa_boot_influx_zt21k_sm10_0500.csv"

curvesfile <- paste(workpath,curvesname,sep="")
outfile <- paste(workpath,"curvepolyfits.csv",sep="")

curvesin <- read.csv(curvesfile)
sc <- curvesin[1:211,]

# regressions
lm.nhemis <- lm(nhemis ~ poly(targ_age,5), data=sc)
summary(lm.nhemis)
plot(sc$targ_age, sc$nhemis, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="NHemis")
lines(sc$targ_age, lm.nhemis$fitted.values, lw=3, col="green")

lm.shemis <- lm(shemis ~ poly(targ_age,5), data=sc)
summary(lm.shemis)
plot(sc$targ_age, sc$shemis, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="SHemis")
lines(sc$targ_age, lm.shemis$fitted.values, lw=3, col="green")

lm.globe <- lm(globe ~ poly(targ_age,5), data=sc)
summary(lm.globe)
plot(sc$targ_age, sc$globe, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="Globe")
lines(sc$targ_age, lm.globe$fitted.values, lw=3, col="green")

lm.ntrop <- lm(ntrop ~ poly(targ_age,5), data=sc)
summary(lm.ntrop)
plot(sc$targ_age, sc$ntrop, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="Northern Tropics")
lines(sc$targ_age, lm.ntrop$fitted.values, lw=3, col="green")

lm.strop <- lm(strop ~ poly(targ_age,5), data=sc)
summary(lm.strop)
plot(sc$targ_age, sc$strop, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="Southern Tropics")
lines(sc$targ_age, lm.strop$fitted.values, lw=3, col="green")

lm.nextr <- lm(nextr ~ poly(targ_age,5), data=sc)
summary(lm.nextr)
plot(sc$targ_age, sc$nextr, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="Northern Extratropics")
lines(sc$targ_age, lm.nextr$fitted.values, lw=3, col="green")

lm.sextr <- lm(sextr ~ poly(targ_age,5), data=sc)
summary(lm.sextr)
plot(sc$targ_age, sc$sextr, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="Southern Extratropics")
lines(sc$targ_age, lm.sextr$fitted.values, lw=3, col="green")

lm.hemis_diff <- lm(hemis_diff ~ poly(targ_age,5), data=sc)
summary(lm.hemis_diff)
plot(sc$targ_age, sc$hemis_diff, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="Hemispheric Difference")
lines(sc$targ_age, lm.hemis_diff$fitted.values, lw=3, col="green")

out <- cbind(scoresin,lm.nhemis$fitted.values,lm.shemis$fitted.values,lm.globe$fitted.values,
    lm.ntrop$fitted.values,lm.strop$fitted.values,lm.nextr$fitted.values,lm.sextr$fitted.values,
    lm.hemis_diff$fitted.values)

write.table(out, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")