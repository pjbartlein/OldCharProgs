# regressions of charcoal curves on PCs

workpath <-"e:\\Projects\\GPWG\\analysis\\output\\v2n\\v2n_curves2\\"
curvesname <- "v2n_all_nsa_boot_influx_zt21k_sm10_0500.csv"

gcmpath <- "e:\\Projects\\Transient\\work01\\"
tmpname <- "tmp2m_land_zave_anm_2000_latzones_lwr.csv"
prename <- "pre_land_zave_anm_2000_latzones_lwr.csv"
pmename <- "p-e_land_zave_anm_2000_latzones_lwr.csv"

curvesfile <- paste(workpath,curvesname,sep="")
outfile <- paste(workpath,"curvegcmpolyfits.csv",sep="")

curvesin <- read.csv(curvesfile)
char <- curvesin[1:211,]
tmp <- read.csv(paste(gcmpath, tmpname, sep=""))
pre <- read.csv(paste(gcmpath, prename, sep=""))
pme <- read.csv(paste(gcmpath, pmename, sep=""))


# regressions
x1 <- tmp$Nhemis
x2 <- pre$Nhemis
x3 <- pme$Nhemis
lm.nhemis <- lm(char$nhemis ~ x1 + x3 +  I(x3^2) )
summary(lm.nhemis)
plot(char$targ_age, char$nhemis, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="NHemis")
lines(char$targ_age, lm.nhemis$fitted.values, lw=3, col="deeppink")

x1 <- tmp$Shemis
x2 <- pre$Shemis
x3 <- pme$Shemis
lm.shemis <- lm(char$shemis ~ x1 + x3 +  I(x3^2) )
plot(char$targ_age, char$shemis, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="SHemis")
lines(char$targ_age, lm.shemis$fitted.values, lw=3, col="deeppink")

x1 <- tmp$Globe
x2 <- pre$Globe
x3 <- pme$Globe
lm.globe <- lm(char$globe ~ x1 + x3 +  I(x3^2) )
summary(lm.globe)
plot(char$targ_age, char$globe, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="Globe")
lines(char$targ_age, lm.globe$fitted.values, lw=3, col="deeppink")

x1 <- tmp$Ntrop
x2 <- pre$Ntrop
x3 <- pme$Ntrop
lm.ntrop <- lm(char$ntrop ~ x1 + x3 +  I(x3^2) )
summary(lm.ntrop)
plot(char$targ_age, char$ntrop, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="Northern Tropics")
lines(char$targ_age, lm.ntrop$fitted.values, lw=3, col="deeppink")

x1 <- tmp$Strop
x2 <- pre$Strop
x3 <- pme$Strop
lm.strop <- lm(char$strop ~ x1 + x3 +  I(x3^2) )
summary(lm.strop)
plot(char$targ_age, char$strop, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="Southern Tropics")
lines(char$targ_age, lm.strop$fitted.values, lw=3, col="deeppink")

x1 <- tmp$Nextr
x2 <- pre$Nextr
x3 <- pme$Nextr
lm.nextr <- lm(char$nextr ~ x1 + x3 +  I(x3^2) )
summary(lm.nextr)
plot(char$targ_age, char$nextr, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="Northern Extratropics")
lines(char$targ_age, lm.nextr$fitted.values, lw=3, col="deeppink")

x1 <- tmp$Sextr
x2 <- pre$Sextr
x3 <- pme$Sextr
lm.sextr <- lm(char$sextr ~ x1 + x3 +  I(x3^2) )
summary(lm.sextr)
plot(char$targ_age, char$sextr, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="Southern Extratropics")
lines(char$targ_age, lm.sextr$fitted.values, lw=3, col="deeppink")


out <- cbind(scoresin,lm.nhemis$fitted.values,lm.shemis$fitted.values,lm.globe$fitted.values,
    lm.ntrop$fitted.values,lm.strop$fitted.values,lm.nextr$fitted.values,lm.sextr$fitted.values)

write.table(out, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")

#plot(x1, lm.ntrop$coefficients["x1"]*x1 + lm.ntrop$coefficients["I(x1^2)"]*x1*x1)
#plot(x2, lm.ntrop$coefficients["x2"]*x2 + lm.ntrop$coefficients["I(x2^2)"]*x2*x2)
#plot(x3, lm.ntrop$coefficients["x3"]*x3 + lm.ntrop$coefficients["I(x3^2)"]*x3*x3)
