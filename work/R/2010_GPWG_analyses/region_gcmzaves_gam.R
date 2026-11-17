# regressions of charcoal curves on PCs
# library(gam)
workpath <-"e:\\Projects\\GPWG\\analysis\\output\\v2n\\v2n_curves2\\"
curvesname <- "v2n_all_nsa_boot_influx_zt21k_sm10_0500.csv"

gcmpath <- "e:\\Projects\\Transient\\work01\\"
tmpname <- "tmp2m_land_zave_anm_2000_latzones_lwr.csv"
prename <- "pre_land_zave_anm_2000_latzones_lwr.csv"
pmename <- "p-e_land_zave_anm_2000_latzones_lwr.csv"

curvesfile <- paste(workpath,curvesname,sep="")
outfile <- paste(workpath,"curvegcmgamfits.csv",sep="")

curvesin <- read.csv(curvesfile)
char <- curvesin[1:211,]
tmp <- read.csv(paste(gcmpath, tmpname, sep=""))
pre <- read.csv(paste(gcmpath, prename, sep=""))
pme <- read.csv(paste(gcmpath, pmename, sep=""))
par(mfcol=c(1,1))

# regressions
x1 <- tmp$Nhemis
x2 <- pre$Nhemis
x3 <- pme$Nhemis
gam.nhemis <- gam(char$nhemis ~ s(x1, df=1) + s(x3, df=3) ) 
summary(gam.nhemis)
plot(char$targ_age, char$nhemis, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="NHemis")
lines(char$targ_age, gam.nhemis$fitted.values, lw=3, col="deeppink")
plot.gam(gam.nhemis, residuals=T, se=T, scale=4)
r <- cor(char$nhemis, gam.nhemis$fitted.values); r2 <- r^2; r; r2
opar <- par(mfcol=c(2,2))
qqnorm(gam.nhemis$residuals); hist(gam.nhemis$residuals);
plot(char$nhemis ~ gam.nhemis$fitted.values); abline(0,1); plot(gam.nhemis$residuals ~ gam.nhemis$fitted.values)
par(opar)

x1 <- tmp$Shemis
x2 <- pre$Shemis
x3 <- pme$Shemis
gam.shemis <- gam(char$shemis ~ s(x1, df=1) + s(x3, df=3) ) 
summary(gam.shemis)
plot(char$targ_age, char$shemis, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="SHemis")
lines(char$targ_age, gam.shemis$fitted.values, lw=3, col="deeppink")
plot.gam(gam.shemis, residuals=T, se=T, scale=4)
r <- cor(char$shemis, gam.shemis$fitted.values); r2 <- r^2; r; r2
opar <- par(mfcol=c(2,2))
qqnorm(gam.shemis$residuals); hist(gam.shemis$residuals);
plot(char$shemis ~ gam.shemis$fitted.values); abline(0,1); plot(gam.shemis$residuals ~ gam.shemis$fitted.values)
par(opar)

x1 <- tmp$Globe
x2 <- pre$Globe
x3 <- pme$Globe
gam.globe <- gam(char$globe ~ s(x1, df=1) + s(x3, df=2) ) 
summary(gam.globe)
plot(char$targ_age, char$globe, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="Globe")
lines(char$targ_age, gam.globe$fitted.values, lw=3, col="deeppink")
plot.gam(gam.globe, residuals=T, se=T, scale=4)
r <- cor(char$globe, gam.globe$fitted.values); r2 <- r^2; r; r2
opar <- par(mfcol=c(2,2))
qqnorm(gam.globe$residuals); hist(gam.globe$residuals);
plot(char$globe ~ gam.globe$fitted.values); abline(0,1); plot(gam.globe$residuals ~ gam.globe$fitted.values)
par(opar)

x1 <- tmp$Ntrop
x2 <- pre$Ntrop
x3 <- pme$Ntrop
gam.ntrop <- gam(char$ntrop ~ s(x1, df=1) +  s(x3, df=4) ) 
summary(gam.ntrop)
plot(char$targ_age, char$ntrop, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="Northern Tropics")
lines(char$targ_age, gam.ntrop$fitted.values, lw=3, col="deeppink")
plot.gam(gam.ntrop, residuals=T, se=T, scale=4)
r <- cor(char$ntrop, gam.ntrop$fitted.values); r2 <- r^2; r; r2
opar <- par(mfcol=c(2,2))
qqnorm(gam.ntrop$residuals); hist(gam.ntrop$residuals);
plot(char$ntrop ~ gam.ntrop$fitted.values); abline(0,1); plot(gam.ntrop$residuals ~ gam.ntrop$fitted.values)
par(opar)

x1 <- tmp$Strop
x2 <- pre$Strop
x3 <- pme$Strop
gam.strop <- gam(char$strop ~ s(x1, df=1) +  s(x3, df=3) ) 
summary(gam.strop)
plot(char$targ_age, char$strop, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="Southern Tropics")
lines(char$targ_age, gam.strop$fitted.values, lw=3, col="deeppink")
plot.gam(gam.strop, residuals=T, se=T, scale=4)
r <- cor(char$strop, gam.strop$fitted.values); r2 <- r^2; r; r2
opar <- par(mfcol=c(2,2))
qqnorm(gam.strop$residuals); hist(gam.strop$residuals);
plot(char$strop ~ gam.strop$fitted.values); abline(0,1); plot(gam.strop$residuals ~ gam.strop$fitted.values)
par(opar)

x1 <- tmp$Nextr
x2 <- pre$Nextr
x3 <- pme$Nextr
gam.nextr <- gam(char$nextr ~ s(x1, df=1) +  s(x3, df=3) ) 
summary(gam.nextr)
plot(char$targ_age, char$nextr, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="Northern Extratropics")
lines(char$targ_age, gam.nextr$fitted.values, lw=3, col="deeppink")
plot.gam(gam.nextr, residuals=T, se=T, scale=4)
r <- cor(char$nhemis, gam.nextr$fitted.values); r2 <- r^2; r; r2
opar <- par(mfcol=c(2,2))
qqnorm(gam.nextr$residuals); hist(gam.nextr$residuals);
plot(char$nextr ~ gam.nextr$fitted.values); abline(0,1); plot(gam.nextr$residuals ~ gam.nextr$fitted.values)
par(opar)

x1 <- tmp$Sextr
x2 <- pre$Sextr
x3 <- pme$Sextr
gam.sextr <- gam(char$sextr ~ s(x1, df=2) +  s(x3, df=3) ) 
summary(gam.sextr)
plot(char$targ_age, char$sextr, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="Southern Extratropics")
lines(char$targ_age, gam.sextr$fitted.values, lw=3, col="deeppink")
plot.gam(gam.sextr, residuals=T, se=T, scale=4)
r <- cor(char$sextr, gam.sextr$fitted.values); r2 <- r^2; r; r2
opar <- par(mfcol=c(2,2))
qqnorm(gam.sextr$residuals); hist(gam.sextr$residuals);
plot(char$sextr ~ gam.sextr$fitted.values); abline(0,1); plot(gam.sextr$residuals ~ gam.sextr$fitted.values)
par(opar)

out <- cbind(scoresin,gam.nhemis$fitted.values,gam.shemis$fitted.values,gam.globe$fitted.values,
    gam.ntrop$fitted.values,gam.strop$fitted.values,gam.nextr$fitted.values,gam.sextr$fitted.values)

write.table(out, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")


new.x <- expand.grid(x1=seq(from=min(x1), to=max(x1), len=40),
                     x3=seq(from=min(x3), to=max(x3), len=40))
gam.sfc <- predict(gam.ntrop, new.x)
persp(gam.sfc, theta=-45)