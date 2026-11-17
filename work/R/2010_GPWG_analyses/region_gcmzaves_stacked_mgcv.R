# regressions of charcoal curves on gcm zonal averages
# library(mgcv)
workpath <-"e:\\Projects\\GPWG\\analysis\\output\\v2n\\v2n_glfiles\\"

climname <- "SIM2b1_all_anmb_0500_lwr.csv"
charname<- "zt_all_0500.csv"
outfile <- paste(workpath, "gamfits_", climname, sep="")
                 
clim <- read.csv(paste(workpath, climname, sep=""))
char <- read.csv(paste(workpath, charname, sep=""))
globe <- data.frame(age=clim$age, zt=char$globe_zt, sefit=char$globe_sefit, seboot=char$globe_seboot,
            tmp=clim$globe_tmp, pre=clim$globe_pre, pme=clim$globe_pme, region="globe")
nhemis <- data.frame(age=clim$age, zt=char$nhemis_zt, sefit=char$nhemis_sefit, seboot=char$nhemis_seboot,
            tmp=clim$nhemis_tmp, pre=clim$nhemis_pre, pme=clim$nhemis_pme, region="nhemis")
shemis <- data.frame(age=clim$age, zt=char$shemis_zt, sefit=char$shemis_sefit, seboot=char$shemis_seboot,
            tmp=clim$shemis_tmp, pre=clim$shemis_pre, pme=clim$shemis_pme, region="shemis")
nextr <- data.frame(age=clim$age, zt=char$nextr_zt, sefit=char$nextr_sefit, seboot=char$nextr_seboot,
            tmp=clim$nextr_tmp, pre=clim$nextr_pre, pme=clim$nextr_pme, region="nextr")
ntrop <- data.frame(age=clim$age, zt=char$ntrop_zt, sefit=char$ntrop_sefit, seboot=char$ntrop_seboot,
            tmp=clim$ntrop_tmp, pre=clim$ntrop_pre, pme=clim$ntrop_pme, region="ntrop")
strop <- data.frame(age=clim$age, zt=char$strop_zt, sefit=char$strop_sefit, seboot=char$strop_seboot,
            tmp=clim$strop_tmp, pre=clim$strop_pre, pme=clim$strop_pme, region="strop")
sextr <- data.frame(age=clim$age, zt=char$sextr_zt, sefit=char$sextr_sefit, seboot=char$sextr_seboot,
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
plot.gam(mgcv.stack, residuals=T, se=T, shade=T, pch=16, cex=0.5, ylim=c(-1,1), ask="no")
#vis.gam(mgcv.stack, plot.type="persp", type="response", color="gray", theta=-60)
gam.check(mgcv.stack)

# anova
mgcv.stack.null <- gam(zt ~ 1, data=stack)
summary(mgcv.stack.null)
anova(mgcv.stack.null, mgcv.stack, test="F")
r <- cor(stack$zt, mgcv.stack$fitted.values); r2 <- r^2; r; r2

# perspective plot
par(mfcol=c(1,1))
new.tmp <- seq(from=-4, to=2, len=40); new.pme <- seq(from=-1.15, to=1.15, len=40)
#new.tmp <- seq(from=min(x1), to=max(x1), len=40); new.pme <- seq(from=min(x3), to=max(x3), len=40)
new.x <- expand.grid(tmp=new.tmp, pme=new.pme)
mgcv.sfc <- predict(mgcv.stack, new.x)
mgcv.sfc <- matrix(mgcv.sfc, 40, 40)
persp(new.tmp, new.pme, mgcv.sfc, theta=-60, phi=20, d=1.5, col="gray",ticktype="detailed",
      zlim=c(-1.25,1), xlab="Temperature", ylab="P-E",
      zlab = "Z-Score of Tranformed Charcoal Influx")

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

#stack$ztp <- ((stack$zt-min(stack$zt))/(max(stack$zt)-min(stack$zt))+.01)^0.33
