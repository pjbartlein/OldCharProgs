# wus influx
# library(gam)

workpath <-"e:\\Projects\\GPWG\\work\\wus\\R\\"
curvesname <- "boot100.csv"
outname <- "boot100_fitted.csv"
outfile <- pass(workpath, outname, sep="")
                 
wusin <- read.csv(paste(workpath, curvesname, sep=""))
summary(wusin)
wus <- wusin[wusin$Year <= 1800,]
#wus <- wusin
summary(wus)

par(mfcol=c(1,1))

# regressions -- 

gam.wus.null <- lm(influx ~ 1, data=wus)
summary(gam.wus.null)

#gam.wus <- gam(influx ~ s(tann, bs="tp", k=2) + s(dai, bs="tp", k=2) + s(tann, dai, bs="tp", k=4) , data=wus )

# tann alone
gam.wus.t <- gam(influx ~ s(tann), data=wus )
summary(gam.wus.t)
plot.gam(gam.wus.t, residuals=T, se=T, pch=16, cex=0.5)
wus.fit <- predict(gam.wus.t, wus)
wusin.fit <- predict(gam.wus.t, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, gam.wus.t$fitted.values); r2 <- r^2; r; r2
anova(gam.wus.null, gam.wus.t,test="F")
AIC(gam.wus.null, gam.wus.t)

# dai alone
gam.wus.d <- gam(influx ~ s(dai), data=wus )
summary(gam.wus.d)
plot.gam(gam.wus.d, residuals=T, se=T, pch=16, cex=0.5)
wus.fit <- predict(gam.wus.d, wus)
wusin.fit <- predict(gam.wus.d, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, gam.wus.d$fitted.values); r2 <- r^2; r; r2
anova(gam.wus.null, gam.wus.d,test="F")
AIC(gam.wus.null, gam.wus.d)

anova(gam.wus.t, gam.wus.d,test="F")
AIC(gam.wus.t, gam.wus.d)

# additive
gam.wus.tda <- gam(influx ~ lo(tann, span=0.6) + lo(dai, span=0.6), data=wus )
summary(gam.wus.tda)
par(mfcol=c(1,1))
plot.gam(gam.wus.tda, residuals=T, se=T, pch=16, cex=0.5)
wus.fit <- predict(gam.wus.tda, wus)
wusin.fit <- predict(gam.wus.tda, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, gam.wus.td$fitted.values); r2 <- r^2; r; r2
anova(gam.wus.null, gam.wus.tda,test="F")
AIC(gam.wus.null, gam.wus.tda)

anova(gam.wus.t, gam.wus.tda,test="F")
AIC(gam.wus.t, gam.wus.tda)

# interaction
gam.wus.td <- gam(influx ~ lo(tann, dai), data=wus )
summary(gam.wus.td)
par(mfcol=c(1,1))
plot.gam(gam.wus.td, residuals=T, se=T, pch=16, cex=0.5)
wus.fit <- predict(gam.wus.td, wus)
wusin.fit <- predict(gam.wus.td, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, gam.wus.td$fitted.values); r2 <- r^2; r; r2
anova(gam.wus.null, gam.wus.td,test="F")
AIC(gam.wus.null, gam.wus.td)

anova(gam.wus.t, gam.wus.td,test="F")
AIC(gam.wus.t, gam.wus.td)

# tann + interaction
gam.wus.ttd <- gam(influx ~ s(tann) +s(tann*dai), data=wus )
summary(gam.wus.ttd)
par(mfcol=c(1,1))
plot.gam(gam.wus.ttd, residuals=T, se=T, pch=16, cex=0.5)
wus.fit <- predict(gam.wus.ttd, wus)
wusin.fit <- predict(gam.wus.ttd, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, gam.wus.ttd$fitted.values); r2 <- r^2; r; r2
anova(gam.wus.null, gam.wus.ttd,test="F")
AIC(gam.wus.null, gam.wus.ttd)

anova(gam.wus.td, gam.wus.ttd,test="F")
AIC(gam.wus.td, gam.wus.ttd)

# dai + interaction
gam.wus.dtd <- gam(influx ~ s(dai) + s(tann*dai), data=wus )
summary(gam.wus.dtd)
par(mfcol=c(1,1))
plot.gam(gam.wus.dtd, residuals=T, se=T, pch=16, cex=0.5)
wus.fit <- predict(gam.wus.dtd, wus)
wusin.fit <- predict(gam.wus.dtd, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, gam.wus.dtd$fitted.values); r2 <- r^2; r; r2
anova(gam.wus.null, gam.wus.dtd,test="F")
AIC(gam.wus.null, gam.wus.dtd)

anova(gam.wus.td, gam.wus.dtd,test="F")
AIC(gam.wus.td, gam.wus.dtd)

anova(gam.wus.ttd, gam.wus.dtd,test="F")
AIC(gam.wus.ttd, gam.wus.dtd)

# full model
gam.wus <- gam(influx ~ s(tann)  + s(dai) + lo(tann, dai, span=0.8), data=wus ) #
summary(gam.wus)
par(mfcol=c(1,1))
plot.gam(gam.wus, residuals=T, se=T, pch=16, cex=0.5)
wus.fit <- predict(gam.wus, wus)
wusin.fit <- predict(gam.wus, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, gam.wus$fitted.values); r2 <- r^2; r; r2
anova(gam.wus.null, gam.wus,test="F")
AIC(gam.wus.null, gam.wus)

#anova(gam.wus.tda, gam.wus,test="F")
#AIC(gam.wus.tda, gam.wus)

# plots
tannmin <- -0.3; tannmax <- 0.15; daimin <- 0.3; daimax <- 0.39; ngrd <- 40
tanninc <- (tannmax-tannmin)/ngrd; daiinc <- (daimax-daimin)/ngrd
new.tann <- seq(from=tannmin, to=tannmax, len=ngrd); new.dai <- seq(from=daimin, to=daimax, len=ngrd)
##new.wus$tmp <- seq(from=min(wus$tmp), to=max(wus$tmp), len=40); new.wus$pme <- seq(from=min(wus$pme), to=max(wus$pme), len=40)
new.x <- expand.grid(tann=new.tann, dai=new.dai)
gam.sfc <- predict(gam.wus, new.x)
gam.sfc <- matrix(gam.sfc, ngrd, ngrd)

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
gam.sfc[mask==0] <- NA

plot(wus$tann, wus$dai, pch=16, col="gray", cex=0.75, xlim=c(-.3,0.15), ylim=c(0.3,0.39),
     xlab="Annual Temperature Anomaly (C)",
     ylab="Drought-Area Index",
      main="Z-Scores of Transformed Charcoal Influx")
      
contour(new.tann, new.dai, gam, add=T,
      col="red", lwd=2, labcex=0.75,
      xlim=c(-.3,.15),
      levels=c(-0.3, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1))



par(mfcol=c(1,1))
image(new.tann, new.dai, gam)
points(wus$tann, wus$dai)

persp(new.tann, new.dai, gam, theta=-35, phi=30, d=1.5, col="gray", ticktype="detailed",
      zlim=c(-0.5,0.5),
      xlab="MAT (C)", ylab="DAI (JJA)",
      zlab = "Z-Score of Transformed Charcoal Influx")

wus.fit <- predict(gam.wus, wus)
wusin.fit <- predict(gam.wus, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")

out.wus.fit <- cbind(wusin$Year,wusin$influx,wusin.fit)
colnames(out.wus.fit) <- c("Year","Influx","Pred")
write.table(out.wus.fit, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")