# regressions of gl$ztcoal curves on PCs
# library(mgcv)
workpath <-"e:\\Projects\\GPWG\\analysis\\output\\v2n\\v2n_glfiles\\"
curvesname <- "z21-0_zgl_0500.csv"
                 
gl <- read.csv(paste(workpath, curvesname, sep=""))
summary(gl)

par(mfcol=c(1,1))
#maxzt <- max(gl$zt); minzt <- min(gl$zt)
#gl$zt2 <-  log10((gl$zt-minzt)/(maxzt-minzt)+0.1)
#plot(gl$zt, gl$zt2)

# regressions

mgcv.null.gl <- gam(zt ~ 1, data=gl)
summary(mgcv.null.gl)

#mgcv.gl <- gam(zt ~ tmp + te(pme, k=4),  data=gl)
#summary(mgcv.gl)

mgcv.gl <- gam(zt ~ te(tmp, bs="tp", k=4) + te(pme, bs="tp", k=8), data=gl )
summary(mgcv.gl)

plot.gam(mgcv.gl, residuals=F, se=T, shade=T, pch=16, cex=0.5, ylim=c(-1,1), xlim=c(-3,3))
vis.gam(mgcv.gl, plot.type="persp", type="response", color="gray", theta=-60)
gam.check(mgcv.gl)
r <- cor(gl$zt, mgcv.gl$fitted.values); r2 <- r^2; r; r2
anova(mgcv.null.gl, mgcv.gl,test="F")

region <- "ntrop"
gcmfile <- paste(workpath, "SIM2b1_", region, "_z.csv", sep="")
outfile <- paste(workpath, "gamfits_", region, ".csv",sep="")
hemis.dat <- read.csv(gcmfile)
zthat <- predict(mgcv.gl, hemis.dat)
plot(hemis.dat$age, zthat, ylim=c(-1,1))

#out <- cbind(age2,nextr,ntrop,strop,sextr,nextr_fit,ntrop_fit,strop_fit,sextr_fit)
#write.table(out, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")

new.tmp <- seq(from=-2, to=2, len=40); new.pme <- seq(from=-2, to=2, len=40)
##new.gl$tmp <- seq(from=min(gl$tmp), to=max(gl$tmp), len=40); new.gl$pme <- seq(from=min(gl$pme), to=max(gl$pme), len=40)
new.x <- expand.grid(tmp=new.tmp, pme=new.pme)
mgcv.sfc <- predict(mgcv.gl, new.x)
mgcv.sfc <- matrix(mgcv.sfc, 40, 40)
par(mfcol=c(1,1))
persp(new.tmp, new.pme, mgcv.sfc, theta=-60, phi=20, d=1.5, col="gray",ticktype="detailed",
      zlim=c(-1,1),
      xlab="Temperature Z-Scores", ylab="P-E Z-Scores",
      zlab = "Z-Score of Tranformed Charcoal Influx")


opar <- par(mfcol=c(2,4))
globe_fit <- predict(mgcv.gl, globe)
plot(char$age, char$globe_zt, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="globe")
points(char$age, globe_fit, pch=16, cex=0.5, col="purple")

nhemis_fit <- predict(mgcv.gl, nhemis)
plot(char$age, char$nhemis_zt, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="nhemis")
points(char$age, nhemis_fit, pch=16, cex=0.5, col="purple")

shemis_fit <- predict(mgcv.gl, shemis)
plot(char$age, char$shemis_zt, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="shemis")
points(char$age, shemis_fit, pch=16, cex=0.5, col="purple")

nextr_fit <- predict(mgcv.gl, nextr)
plot(char$age, char$nextr_zt, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="nextr")
points(char$age, nextr_fit, pch=16, cex=0.5, col="purple")

ntrop_fit <- predict(mgcv.gl, ntrop)
plot(char$age, char$ntrop_zt, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="ntrop")
points(char$age, ntrop_fit, pch=16, cex=0.5, col="purple")

strop_fit <- predict(mgcv.gl, strop)
plot(char$age, char$strop_zt, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="strop")
points(char$age, strop_fit, pch=16, cex=0.5, col="purple")

sextr_fit <- predict(mgcv.gl, sextr)
plot(char$age, char$sextr_zt, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="sextr")
points(char$age, sextr_fit, pch=16, cex=0.5, col="purple")
par(opar)
