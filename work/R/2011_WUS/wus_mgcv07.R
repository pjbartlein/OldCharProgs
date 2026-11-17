# wus influx
# library(mgcv)

workpath <-"e:\\Projects\\GPWG\\work\\wus\\R\\"
curvesname <- "boot100.csv"
outname <- "boot100_fitted.csv"
outfile <- paste(workpath, outname, sep="")
                 
wusin <- read.csv(paste(workpath, curvesname, sep=""))
summary(wusin)
wus <- wusin[wusin$Year <= 1800,]
#wus <- wusin
summary(wus)

par(mfcol=c(1,1))

# regressions -- 

mgcv.wus.null <- lm(influx ~ 1, data=wus)
summary(mgcv.wus.null)

#mgcv.wus <- gam(influx ~ te(tann, bs="tp", k=2) + te(dai, bs="tp", k=2) + te(tann, dai, bs="tp", k=4) , data=wus )

# tann alone
mgcv.wus.t <- gam(influx ~ s(tann, bs="tp", fx=T, k=2), data=wus )
summary(mgcv.wus.t)
plot.gam(mgcv.wus.t, residuals=T, se=T, shade=T, pch=16, cex=0.5)
gam.check(mgcv.wus.t)
wus.fit <- predict(mgcv.wus.t, wus)
wusin.fit <- predict(mgcv.wus.t, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, mgcv.wus.t$fitted.values); r2 <- r^2; r; r2
anova(mgcv.wus.null, mgcv.wus.t,test="F")
AIC(mgcv.wus.null, mgcv.wus.t)

# dai alone
mgcv.wus.d <- gam(influx ~ s(dai, bs="tp", fx=T, k=4), data=wus )
summary(mgcv.wus.d)
plot.gam(mgcv.wus.d, residuals=T, se=T, shade=T, pch=16, cex=0.5)
gam.check(mgcv.wus.d)
wus.fit <- predict(mgcv.wus.d, wus)
wusin.fit <- predict(mgcv.wus.d, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, mgcv.wus.d$fitted.values); r2 <- r^2; r; r2
anova(mgcv.wus.null, mgcv.wus.d,test="F")
AIC(mgcv.wus.null, mgcv.wus.d)

anova(mgcv.wus.t, mgcv.wus.d,test="F")
AIC(mgcv.wus.t, mgcv.wus.d)

# additive
mgcv.wus.tda <- gam(influx ~ te(tann, bs="tp", fx=T, k=2) + te(dai, bs="tp", fx=T, k=3), data=wus )
summary(mgcv.wus.tda)
par(mfcol=c(1,1))
plot.gam(mgcv.wus.tda, residuals=T, se=T, shade=T, pch=16, cex=0.5)
vis.gam(mgcv.wus, plot.type="persp", type="response", color="gray", theta=-60)
gam.check(mgcv.wus.tda)
wus.fit <- predict(mgcv.wus.tda, wus)
wusin.fit <- predict(mgcv.wus.tda, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, mgcv.wus.td$fitted.values); r2 <- r^2; r; r2
anova(mgcv.wus.null, mgcv.wus.tda,test="F")
AIC(mgcv.wus.null, mgcv.wus.tda)

anova(mgcv.wus.t, mgcv.wus.tda,test="F")
AIC(mgcv.wus.t, mgcv.wus.tda)

# interaction
mgcv.wus.td <- gam(influx ~ s(tann, dai, bs="tp", fx=T, k=4), data=wus )
summary(mgcv.wus.td)
par(mfcol=c(1,1))
plot.gam(mgcv.wus.td, residuals=T, se=T, shade=T, pch=16, cex=0.5)
vis.gam(mgcv.wus, plot.type="persp", type="response", color="gray", theta=-60)
gam.check(mgcv.wus.td)
wus.fit <- predict(mgcv.wus.td, wus)
wusin.fit <- predict(mgcv.wus.td, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, mgcv.wus.td$fitted.values); r2 <- r^2; r; r2
anova(mgcv.wus.null, mgcv.wus.td,test="F")
AIC(mgcv.wus.null, mgcv.wus.td)

anova(mgcv.wus.t, mgcv.wus.td,test="F")
AIC(mgcv.wus.t, mgcv.wus.td)

# tann + interaction
mgcv.wus.ttd <- gam(influx ~ te(tann, bs="tp", fx=T, k=2) +te(tann, dai, bs="tp", fx=F, k=4), data=wus )
summary(mgcv.wus.ttd)
par(mfcol=c(1,1))
plot.gam(mgcv.wus.ttd, residuals=T, se=T, shade=T, pch=16, cex=0.5)
vis.gam(mgcv.wus, plot.type="persp", type="response", color="gray", theta=-60)
gam.check(mgcv.wus.ttd)
wus.fit <- predict(mgcv.wus.ttd, wus)
wusin.fit <- predict(mgcv.wus.ttd, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, mgcv.wus.ttd$fitted.values); r2 <- r^2; r; r2
anova(mgcv.wus.null, mgcv.wus.ttd,test="F")
AIC(mgcv.wus.null, mgcv.wus.ttd)

anova(mgcv.wus.td, mgcv.wus.ttd,test="F")
AIC(mgcv.wus.td, mgcv.wus.ttd)

# dai + interaction
mgcv.wus.dtd <- gam(influx ~ te(dai, bs="tp", fx=F, k=2) + te(tann, dai, bs="tp", fx=F, k=4), data=wus )
summary(mgcv.wus.dtd)
par(mfcol=c(1,1))
plot.gam(mgcv.wus.dtd, residuals=T, se=T, shade=T, pch=16, cex=0.5)
gam.check(mgcv.wus.dtd)
wus.fit <- predict(mgcv.wus.dtd, wus)
wusin.fit <- predict(mgcv.wus.dtd, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, mgcv.wus.dtd$fitted.values); r2 <- r^2; r; r2
anova(mgcv.wus.null, mgcv.wus.dtd,test="F")
AIC(mgcv.wus.null, mgcv.wus.dtd)

anova(mgcv.wus.td, mgcv.wus.dtd,test="F")
AIC(mgcv.wus.td, mgcv.wus.dtd)

anova(mgcv.wus.ttd, mgcv.wus.dtd,test="F")
AIC(mgcv.wus.ttd, mgcv.wus.dtd)

# full model
mgcv.wus <- gam(influx ~ te(tann, bs="tp", fx=T, k=2)  + te(dai, bs="tp", fx=T, k=3) + te(tann, dai, bs="tp", fx=F, k=2), data=wus ) #
summary(mgcv.wus)
par(mfcol=c(1,1))
plot.gam(mgcv.wus, residuals=T, se=F, shade=T, pch=16, cex=0.5)
vis.gam(mgcv.wus, plot.type="persp", type="response", color="gray", theta=-60)
gam.check(mgcv.wus)
wus.fit <- predict(mgcv.wus, wus)
wusin.fit <- predict(mgcv.wus, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, mgcv.wus$fitted.values); r2 <- r^2; r; r2
anova(mgcv.wus.null, mgcv.wus,test="F")
AIC(mgcv.wus.null, mgcv.wus)

#anova(mgcv.wus.tda, mgcv.wus,test="F")
#AIC(mgcv.wus.tda, mgcv.wus)

# plots
tannmin <- -0.3; tannmax <- 0.15; daimin <- 0.3; daimax <- 0.39; ngrd <- 40
tanninc <- (tannmax-tannmin)/ngrd; daiinc <- (daimax-daimin)/ngrd
new.tann <- seq(from=tannmin, to=tannmax, len=ngrd); new.dai <- seq(from=daimin, to=daimax, len=ngrd)
new.x <- expand.grid(tann=new.tann, dai=new.dai)
mgcv <- predict(mgcv.wus, new.x)
mgcv <- matrix(mgcv, ngrd, ngrd)

mask <- table(cut(wus$tann,breaks=c(seq(from=tannmin-tanninc/2, to=tannmax+tanninc/2, len=ngrd+1))),
              cut(wus$dai,breaks=c(seq(from=daimin-daiinc/2, to=daimax+daiinc/2, len=ngrd+1))))
mask <- as.matrix(mask)
mask[mask > 0] <- 1
#image(mask)

# spread the masked cells
for (k in seq(1,5,1)) {
      mask2 <- mask
      for (i in seq(1,ngrd-1,1)) {
            for (j in seq(1,ngrd-1,1)) {
                 if (mask[i,j]==1) { mask2[i,j+1] <- 1; mask2[i+1,j] <- 1}
      }}
      for (i in seq(2,ngrd,1)) {
            for (j in seq(2,ngrd-1,1)) {
                 if (mask[i,j]==1) { mask2[i,j-1] <- 1; mask2[i-1,j] <- 1}
      }}
      mask <- mask2      
#image(mask)
}
mgcv[mask==0] <- NA

plot(wus$tann, wus$dai, pch=16, col="gray", cex=0.75, xlim=c(-.3,0.15), ylim=c(0.3,0.39),
     xlab="Annual Temperature Anomaly (C)",
     ylab="Drought-Area Index",
      main="Z-Scores of Transformed Charcoal Influx")
      
contour(new.tann, new.dai, mgcv, add=T,
      col="red", lwd=2, labcex=0.75,
      xlim=c(-.3,.15),
      levels=c(-0.3, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1))



par(mfcol=c(1,1))
image(new.tann, new.dai, mgcv)
points(wus$tann, wus$dai)

persp(new.tann, new.dai, mgcv, theta=-60, phi=10, d=1.5, col="gray", ticktype="detailed",
      zlim=c(-0.5,0.5),
      xlab="MAT (C)", ylab="DAI (JJA)",
      zlab = "Z-Score of Transformed Charcoal Influx")

wus.fit <- predict(mgcv.wus, wus)
wusin.fit <- predict(mgcv.wus, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")

out.wus.fit <- cbind(wusin$Year,wusin$influx,wusin.fit)
colnames(out.wus.fit) <- c("Year","Influx","Pred")
write.table(out.wus.fit, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")