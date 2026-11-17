# regressions of charcoal curves on gcm zonal averages
# library(mgcv)
workpath <-"e:\\Projects\\GPWG\\analysis\\output\\v2n\\v2n_glfiles\\"

climname <- "SIM2b1_all_anmb_0500_lwr.csv"
charname<- "zt_all_0500.csv"
outfile <- paste(workpath, "finalfits_", climname, sep="")
outfile2 <- paste(workpath, "finalfits_1d_", climname, sep="")
                 
clim <- read.csv(paste(workpath, climname, sep=""))
char <- read.csv(paste(workpath, charname, sep=""))
globe <- data.frame(age=clim$age, zt=char$globe_zt, 
            tmp=clim$globe_tmp, pre=clim$globe_pre, pme=clim$globe_pme, region="globe")
nhemis <- data.frame(age=clim$age, zt=char$nhemis_zt, 
            tmp=clim$nhemis_tmp, pre=clim$nhemis_pre, pme=clim$nhemis_pme, region="nhemis")
shemis <- data.frame(age=clim$age, zt=char$shemis_zt, 
            tmp=clim$shemis_tmp, pre=clim$shemis_pre, pme=clim$shemis_pme, region="shemis")
nextr <- data.frame(age=clim$age, zt=char$nextr_zt, 
            tmp=clim$nextr_tmp, pre=clim$nextr_pre, pme=clim$nextr_pme, region="nextr")
ntrop <- data.frame(age=clim$age, zt=char$ntrop_zt, 
            tmp=clim$ntrop_tmp, pre=clim$ntrop_pre, pme=clim$ntrop_pme, region="ntrop")
strop <- data.frame(age=clim$age, zt=char$strop_zt, 
            tmp=clim$strop_tmp, pre=clim$strop_pre, pme=clim$strop_pme, region="strop")
sextr <- data.frame(age=clim$age, zt=char$sextr_zt,
            tmp=clim$sextr_tmp, pre=clim$sextr_pre, pme=clim$sextr_pme, region="sextr")

nextr.in <- nextr[nextr$age <= 21000,]
ntrop.in <- ntrop[ntrop$age <= 21000,]
strop.in <- strop[strop$age <= 21000,]
sextr.in <- sextr[sextr$age <= 21000,]

stack <-rbind(nextr.in, ntrop.in, strop.in, sextr.in)
plot(stack$age, stack$pme, col=stack$region)
plot(stack$age, stack$tmp, col=stack$region)

# regressions
# final model
mgcv.stack <- gam(zt ~te(tmp, bs="tp", k=3.5) + te(pme, bs="tp", k=8), data=stack)
summary(mgcv.stack)

# diagnostic plots
par(mfcol=c(1,1))
plot.gam(mgcv.stack, residuals=F, se=T, shade=T, pch=16, cex=0.5, ylim=c(-1,1), ask="no")
#vis.gam(mgcv.stack, plot.type="persp", type="response", color="gray", theta=-60)
gam.check(mgcv.stack)

# anova
mgcv.stack.null <- gam(zt ~ 1, data=stack)
summary(mgcv.stack.null)
anova(mgcv.stack.null, mgcv.stack, test="F")
r <- cor(stack$zt, mgcv.stack$fitted.values); r2 <- r^2; r; r2

# perspective plot
par(mfcol=c(1,1))
new.tmp <- seq(from=-3, to=3, len=40); new.pme <- seq(from=-1.15, to=1.15, len=40)
#new.tmp <- seq(from=min(x1), to=max(x1), len=40); new.pme <- seq(from=min(x3), to=max(x3), len=40)
new.x <- expand.grid(tmp=new.tmp, pme=new.pme)
mgcv.sfc <- predict(mgcv.stack, new.x)
mgcv.sfc <- matrix(mgcv.sfc, 40, 40)
persp(new.tmp, new.pme, mgcv.sfc, theta=-65, phi=30, d=1.5, col="gray",ticktype="detailed", #-60 20
      zlim=c(-1.25,1), xlab="MAT (Z-Scores)", ylab="Annual P-E (Z-Scores)",
      zlab = "Z-Scores of Transformed Charcoal Influx")

opar <- par(mfcol=c(2,4))
globe_fit <- predict(mgcv.stack, globe)
plot(char$age, char$globe_zt, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="globe")
points(char$age, globe_fit, pch=16, cex=0.5, col="purple")

nhemis_fit <- predict(mgcv.stack, nhemis)
plot(char$age, char$nhemis_zt, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="nhemis")
points(char$age, nhemis_fit, pch=16, cex=0.5, col="purple")

shemis_fit <- predict(mgcv.stack, shemis)
plot(char$age, char$shemis_zt, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="shemis")
points(char$age, shemis_fit, pch=16, cex=0.5, col="purple")

nextr_fit <- predict(mgcv.stack, nextr)
plot(char$age, char$nextr_zt, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="nextr")
points(char$age, nextr_fit, pch=16, cex=0.5, col="purple")

ntrop_fit <- predict(mgcv.stack, ntrop)
plot(char$age, char$ntrop_zt, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="ntrop")
points(char$age, ntrop_fit, pch=16, cex=0.5, col="purple")

strop_fit <- predict(mgcv.stack, strop)
plot(char$age, char$strop_zt, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="strop")
points(char$age, strop_fit, pch=16, cex=0.5, col="purple")

sextr_fit <- predict(mgcv.stack, sextr)
plot(char$age, char$sextr_zt, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="sextr")
points(char$age, sextr_fit, pch=16, cex=0.5, col="purple")
par(opar)

out <- cbind(globe$age, globe$zt,nhemis$zt,shemis$zt,nextr$zt,ntrop$zt,strop$zt,sextr$zt,
             globe_fit,nhemis_fit,shemis_fit,nextr_fit,ntrop_fit,strop_fit,sextr_fit)
colnames(out) <- c("age","globe_zt","nhemis_zt","shemis_zt","nextr_zt","ntrop_zt","strop_zt","sextr_zt",
             "globe_fit","nhemis_fit","shemis_fit","nextr_fit","ntrop_fit","strop_fit","sextr_fit")
write.table(out, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")

new.tmp <- seq(from=-3, to=3, len=40); new.pme <- seq(from=-1.0, to=1.0, len=40)
new.pme2 <- 0.0
new.x <- expand.grid(tmp=new.tmp, pme=new.pme2)
mgcv.fit.tmp <- predict(mgcv.stack, new.x, se.fit="TRUE")
new.tmp2 <- 1.0
new.x <- expand.grid(tmp=new.tmp2, pme=new.pme)
mgcv.fit.pme <- predict(mgcv.stack, new.x, se.fit="TRUE")

mgcv.stack.1dse <- cbind(new.tmp, mgcv.fit.tmp$fit, mgcv.fit.tmp$fit+2*mgcv.fit.tmp$se.fit, mgcv.fit.tmp$fit-2*mgcv.fit.tmp$se.fit,
                       new.pme, mgcv.fit.pme$fit, mgcv.fit.pme$fit+2*mgcv.fit.pme$se.fit, mgcv.fit.pme$fit-2*mgcv.fit.pme$se.fit) 
colnames(mgcv.stack.1dse) <- c("tmp","tmpfit","tmpfitlcl","tmpfitucl","pme","pmefit","pmelcl","pmeucl")
write.table(mgcv.stack.1dse, file=outfile2, row.names=FALSE, col.names=TRUE, sep=",")



mgcv.nextr <- gam(zt ~te(tmp, bs="tp", k=3.5) + te(pme, bs="tp", k=8), data=nextr)
summary(mgcv.nextr)

# diagnostic plots
par(mfcol=c(1,1))
plot.gam(mgcv.nextr, residuals=F, se=T, shade=T, pch=16, cex=0.5, ylim=c(-1,1), ask="no")
#vis.gam(mgcv.stack, plot.type="persp", type="response", color="gray", theta=-60)
gam.check(mgcv.nextr)

# anova
mgcv.nextr.null <- gam(zt ~ 1, data=nextr)
summary(mgcv.nextr.null)
anova(mgcv.nextr.null, mgcv.nextr, test="F")
r <- cor(nextr$zt, mgcv.nextr$fitted.values); r2 <- r^2; r; r2

# perspective plot
par(mfcol=c(1,1))
new.tmp <- seq(from=-4, to=4, len=40); new.pme <- seq(from=-2, to=2.0, len=40)
#new.tmp <- seq(from=min(x1), to=max(x1), len=40); new.pme <- seq(from=min(x3), to=max(x3), len=40)
new.x <- expand.grid(tmp=new.tmp, pme=new.pme)
mgcv.sfc <- predict(mgcv.nextr, new.x)
mgcv.sfc <- matrix(mgcv.sfc, 40, 40)
persp(new.tmp, new.pme, mgcv.sfc, theta=-65, phi=30, d=1.5, col="gray",ticktype="detailed", #-60 20
      zlim=c(-2,2), xlab="MAT Anomaly", ylab="Annual P-E Anomaly",
      zlab = "Z-Scores of N. Extratropics Transformed Charcoal Influx")

outfile2 <- paste(workpath, "nextr_fits_1d_", climname, sep="")

new.tmp <- seq(from=-4, to=4, len=81); new.pme <- seq(from=-2.0, to=2.0, len=81)
new.pme2 <- 0.0
new.x <- expand.grid(tmp=new.tmp, pme=new.pme2)
mgcv.fit.tmp <- predict(mgcv.nextr, new.x, se.fit="TRUE")
new.tmp2 <- 1.0
new.x <- expand.grid(tmp=new.tmp2, pme=new.pme)
mgcv.fit.pme <- predict(mgcv.nextr, new.x, se.fit="TRUE")

mgcv.nextr.1dse <- cbind(new.tmp, mgcv.fit.tmp$fit, mgcv.fit.tmp$fit+2*mgcv.fit.tmp$se.fit, mgcv.fit.tmp$fit-2*mgcv.fit.tmp$se.fit,
                       new.pme, mgcv.fit.pme$fit, mgcv.fit.pme$fit+2*mgcv.fit.pme$se.fit, mgcv.fit.pme$fit-2*mgcv.fit.pme$se.fit) 
colnames(mgcv.nextr.1dse) <- c("tmp","tmpfit","tmpfitlcl","tmpfitucl","pme","pmefit","pmelcl","pmeucl")
write.table(mgcv.nextr.1dse, file=outfile2, row.names=FALSE, col.names=TRUE, sep=",")