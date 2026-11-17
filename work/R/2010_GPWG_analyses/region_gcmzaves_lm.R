# regressions of charcoal curves on PCs
workpath <-"e:\\Projects\\GPWG\\analysis\\output\\v2n\\v2n_curves2\\"
curvesname <- "v2n_all_nsa_boot_influx_zt21k_sm10_0500.csv"

gcmpath <- "e:\\Projects\\Transient\\work01\\"
tmpname <- "tmp2m_land_zave_anm_2000_latzones_lwr.csv"
prename <- "pre_land_zave_anm_2000_latzones_lwr.csv"
pmename <- "p-e_land_zave_anm_2000_latzones_lwr.csv"

curvesfile <- paste(workpath,curvesname,sep="")
outfile <- paste(workpath,"curvetmpfits.csv",sep="")

curvesin <- read.csv(curvesfile)
char <- curvesin[1:211,]
tmp <- read.csv(paste(gcmpath, tmpname, sep=""))
pre <- read.csv(paste(gcmpath, prename, sep=""))
pme <- read.csv(paste(gcmpath, pmename, sep=""))


# regressions
lm.nhemis <- lm(char$nhemis ~ tmp$Nhemis )
summary(lm.nhemis)
plot(char$targ_age, char$nhemis, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="NHemis")
lines(char$targ_age, lm.nhemis$fitted.values, lw=3, col="orange")

lm.shemis <- lm(char$shemis ~ tmp$Shemis )
summary(lm.shemis)
plot(char$targ_age, char$shemis, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="SHemis")
lines(char$targ_age, lm.shemis$fitted.values, lw=3, col="orange")

lm.globe <- lm(char$globe ~ tmp$Globe )
summary(lm.globe)
plot(char$targ_age, char$globe, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="Globe")
lines(char$targ_age, lm.globe$fitted.values, lw=3, col="orange")

lm.ntrop <- lm(char$ntrop ~ tmp$Ntrop )
summary(lm.ntrop)
plot(char$targ_age, char$ntrop, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="Northern Tropics")
lines(char$targ_age, lm.ntrop$fitted.values, lw=3, col="orange")

lm.strop <- lm(char$strop ~ tmp$Strop )
summary(lm.strop)
plot(char$targ_age, char$strop, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="Southern Tropics")
lines(char$targ_age, lm.strop$fitted.values, lw=3, col="orange")

lm.nextr <- lm(char$nextr ~ tmp$Nextr )
summary(lm.nextr)
plot(char$targ_age, char$nextr, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="Northern Extratropics")
lines(char$targ_age, lm.nextr$fitted.values, lw=3, col="orange")

lm.sextr <- lm(char$sextr ~ tmp$Sextr )
summary(lm.sextr)
plot(char$targ_age, char$sextr, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="Southern Extratropics")
lines(char$targ_age, lm.sextr$fitted.values, lw=3, col="orange")


out <- cbind(scoresin,lm.nhemis$fitted.values,lm.shemis$fitted.values,lm.globe$fitted.values,
    lm.ntrop$fitted.values,lm.strop$fitted.values,lm.nextr$fitted.values,lm.sextr$fitted.values)

write.table(out, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")

AIC(lm.nhemis,gam.nhemis)
anova(lm.nhemis,gam.nhemis)
AIC(lm.shemis,gam.shemis)
anova(lm.shemis,gam.shemis)
AIC(lm.globe,gam.globe)
anova(lm.globe,gam.globe)
AIC(lm.ntrop,gam.ntrop)
anova(lm.ntrop,gam.ntrop)
AIC(lm.strop,gam.strop)
anova(lm.strop,gam.strop)
AIC(lm.nextr,gam.nextr)
anova(lm.nextr,gam.nextr)
AIC(lm.nextr,gam.sextr)
anova(lm.sextr,gam.sextr)