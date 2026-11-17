# regressions of global points
# library(mgcv)

workpath <-"e:\\Projects\\GPWG\\analysis\\output\\v2n\\v2n_glfiles\\"
curvesname <- "z21-0_zgl_0500.csv"

gcmname <- "SIM2b1_all_z_0500_lwr.csv"
outname <- "z21-0_zgl_0500_fitted.csv"
outfile <- paste(workpath, outname, sep="")
                 
gl <- read.csv(paste(workpath, curvesname, sep=""))
summary(gl)

par(mfcol=c(1,1))

# regressions

mgcv.null.gl <- gam(zt ~ 1, data=gl)
summary(mgcv.null.gl)

mgcv.gl <- gam(zt ~ te(tmp, bs="tp", k=4) + te(pme, bs="tp", k=8), data=gl )
summary(mgcv.gl)

plot.gam(mgcv.gl, residuals=F, se=T, shade=T, pch=16, cex=0.5, ylim=c(-1,1), xlim=c(-3,3))
#vis.gam(mgcv.gl, plot.type="persp", type="response", color="gray", theta=-60)
gam.check(mgcv.gl)
r <- cor(gl$zt, mgcv.gl$fitted.values); r2 <- r^2; r; r2
anova(mgcv.null.gl, mgcv.gl,test="F")

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

gcm <- read.csv(paste(workpath, gcmname, sep=""))
summary(gcm)
new.gl <- data.frame(cbind(tmp=gcm$globe_tmp, pme=gcm$globe_pme))

gl.fit <- predict(mgcv.gl, new.gl)

out.gl.fit <- cbind(gcm$age,gl.fit)
write.table(out.gl.fit, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")