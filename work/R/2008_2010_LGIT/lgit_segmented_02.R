# piecewise regression fits
library(segmented)
library(splines)

workpath <-"e:\\Projects\\GPWG\\analysis\\output\\v2n\\v2n_curves2\\"
#curvesname <- "v2n_na_presample_influx_zt21k.csv"
curvesname <- "v2n_eur_presample_influx_zt21k.csv"
#curvesname <- "v2n_ssa_presample_influx_zt21k.csv"
#curvesname <- "v2n_smainaustl_presample_influx_zt21k.csv"
#curvesname <- "v2n_saustlnz_presample_influx_zt21k.csv"

infile <- paste(workpath,curvesname,sep="")
#outfile <- paste(workpath,"v2n_na_presample_influx_zt21k_seg.csv",sep="")
outfile <- paste(workpath,"v2n_eur_presample_influx_zt21k_seg2.csv",sep="")
#outfile <- paste(workpath,"v2n_ssa_presample_influx_zt21k_seg.csv",sep="")
#outfile <- paste(workpath,"v2n_smainaustl_presample_influx_zt21k_seg.csv",sep="")
#outfile <- paste(workpath,"v2n_saustlnz_presample_influx_zt21k_seg.csv",sep="")

lgitin <- read.csv(infile)
lgit <- lgitin[lgitin$age >= 10000.0 & lgitin$age <= 15000.0,]
age2 <- seq(10000,15000, by=50)

# basic trend
lm_zt <- lm(zt ~ age, data=lgit)
summary(lm_zt)
plot(lgit$zt ~ lgit$age, pch=16, cex=0.5, col="gray", xlim=c(15000,10000), xlab="Age", ylab="zt")
points(predict(lm_zt, data.frame(age=age2)) ~ age2, pch=16, cex=0.5, col="red")

# segmented trend
start_breaks <- c(13000, 12000) # c(14000, 11500) # c(14200, 12600, 11800) # 
seg_zt <- segmented(lm_zt, seg.Z = ~ age, psi=list(age=start_breaks), 
	control=seg.control(toll=1e-4, display=T, stop.if.error=F, it.max=500, maxit.glm=20, h=0.1))
summary(seg_zt)
slope(seg_zt)

# plot segmented trend
plot(lgit$zt ~ lgit$age, pch=16, cex=0.5, col="gray", xlim=c(15000,10000), xlab="Age", ylab="zt")
plot(seg_zt, add=T, col="purple", show.gap=F, lwd=3)
#points(seg_zt$fitted.values ~ lgit$age, col="purple", pch=16, cex=0.5)
lines(seg_zt, col="purple", lwd=2)

# piecewise regression via splines package -- segmented breakpoints
knots.seg <- seg_zt$psi[,2]
lmbs_zt <- lm(zt ~ bs(age, knots=knots.seg, degree=1), data=lgit)
summary(lmbs_zt)
plot(lgit$zt ~ lgit$age, pch=16, cex=0.5, col="gray", xlim=c(15000,10000), xlab="Age", ylab="zt")
pred_lmbs_zt <- predict(lmbs_zt, data.frame(age=age2))
points(pred_lmbs_zt ~ age2, pch=16, cex=0.5, col="blue")

# piecewise regression via splines package -- lgit breakpoints
knots.lgit <- c(14200, 12600) # c( 14200, 12600, 11800) # c(12900, 11700) 
lmbs_zt2 <- lm(zt ~ bs(age, knots=knots.lgit, degree=1), data=lgit)
summary(lmbs_zt2)
plot(lgit$zt ~ lgit$age, pch=16, cex=0.5, col="gray", xlim=c(15000,10000), xlab="Age", ylab="zt")
pred_lmbs_zt2 <- predict(lmbs_zt2, data.frame(age=age2))
points(pred_lmbs_zt2  ~ age2, pch=16, cex=0.5, col="green")

AIC(lm_zt,seg_zt,lmbs_zt,lmbs_zt2)

AIC#output
seg_out <- data.frame(cbind(age2,pred_lmbs_zt,pred_lmbs_zt2))
names(seg_out) <- c("year", "lmbs_zt", "lmbs_zt2")
write.table(seg_out, outfile, row.names=FALSE, sep=",")
