# regressions of gl$ztcoal curves on PCs
# 
workpath <-"e:\\Projects\\GPWG\\analysis\\output\\v2n\\v2n_glfiles\\"
curvesname <- "z21-0_zgl_0500.csv"
                 
gl <- read.csv(paste(workpath, curvesname, sep=""))
summary(gl)

par(mfcol=c(1,1))
#maxzt <- max(gl$zt); minzt <- min(gl$zt)
#gl$zt2 <-  log10((gl$zt-minzt)/(maxzt-minzt)+0.1)
#plot(gl$zt, gl$zt2)

# regressions
#x1 + x2 + x3 + I(x1*x2) + I(x1*x3) + I(x2*x3) + I(x1^2) + I(x2^2) + I(x3^2)

lm.gl <- lm(zt ~ tmp + pme + I(tmp*pme) + I(tmp^2) + I(pme^2),  data=gl)
summary(lm.gl)

anova(lm.gl, mgcv.gl,test="F")

#plot(lm.gl)


region <- "strop"
gcmfile <- paste(workpath, "SIM2b1_", region, "_zaves.csv", sep="")
outfile <- paste(workpath, "gamfits_", region, ".csv",sep="")
hemis.dat <- read.csv(gcmfile)
zthat <- predict(lm.gl, hemis.dat)
plot(hemis.dat$age, zthat, ylim=c(-1,1))

#out <- cbind(age2,nextr,ntrop,strop,sextr,nextr_fit,ntrop_fit,strop_fit,sextr_fit)
#write.table(out, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")

new.tmp <- seq(from=-2, to=2, len=40); new.pme <- seq(from=-2, to=2, len=40)
##new.gl$tmp <- seq(from=min(gl$tmp), to=max(gl$tmp), len=40); new.gl$pme <- seq(from=min(gl$pme), to=max(gl$pme), len=40)
new.x <- expand.grid(tmp=new.tmp, pme=new.pme)
mgcv.sfc <- predict(lm.gl, new.x)
mgcv.sfc <- matrix(mgcv.sfc, 40, 40)
par(mfcol=c(1,1))
persp(new.tmp, new.pme, mgcv.sfc, theta=-60, phi=20, d=1.5, col="gray",ticktype="detailed",
      zlim=c(-1,1),
      xlab="Temperature X-Score", ylab="P-E Z-Score",
      zlab = "Z-Score of Tranformed Charcoal Influx")
