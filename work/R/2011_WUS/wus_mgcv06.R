# wus influx
# library(mgcv)

workpath <-"e:\\Projects\\GPWG\\work\\wus\\R\\"
curvesname <- "boot100.csv"
outname <- "boot100p_fitted.csv"
outfile <- paste(workpath, outname, sep="")
                 
wusin <- read.csv(paste(workpath, curvesname, sep=""))
summary(wusin)
wus <- wusin[wusin$Year <= 1800,]
#wus <- wusin
summary(wus)

par(mfcol=c(1,1))

# regressions

mgcv.wus.null <- gam(influx ~ 1, data=wus)
summary(mgcv.wus.null)

#mgcv.wus <- gam(influx ~ te(tann, bs="tp", k=2) + te(pdsi, bs="tp", k=2) + te(tann, pdsi, bs="tp", k=4) , data=wus )

# tann alone
mgcv.wus.t <- gam(influx ~ te(tann, bs="tp", fx=T, k=4), data=wus )
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

# pdsi alone
mgcv.wus.p <- gam(influx ~ te(pdsi, bs="tp", fx=T, k=4), data=wus )
summary(mgcv.wus.p)
plot.gam(mgcv.wus.p, residuals=T, se=T, shade=T, pch=16, cex=0.5)
gam.check(mgcv.wus.p)
wus.fit <- predict(mgcv.wus.p, wus)
wusin.fit <- predict(mgcv.wus.p, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, mgcv.wus.p$fitted.values); r2 <- r^2; r; r2
anova(mgcv.wus.null, mgcv.wus.p,test="F")
AIC(mgcv.wus.null, mgcv.wus.p)

anova(mgcv.wus.t, mgcv.wus.p,test="F")
AIC(mgcv.wus.t, mgcv.wus.p)

# additive
mgcv.wus.tpa<- gam(influx ~ te(tann, bs="tp", fx=T, k=4) + te(pdsi, bs="tp", fx=T, k=4), data=wus )
summary(mgcv.wus.tpa)
par(mfcol=c(1,1))
plot.gam(mgcv.wus.tpa, residuals=T, se=T, shade=T, pch=16, cex=0.5)
vis.gam(mgcv.wus, plot.type="persp", type="response", color="gray", theta=-60)
gam.check(mgcv.wus.tpa)
wus.fit <- predict(mgcv.wus.tpa, wus)
wusin.fit <- predict(mgcv.wus.tpa, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, mgcv.wus.tpa$fitted.values); r2 <- r^2; r; r2
anova(mgcv.wus.null, mgcv.wus.tpa,test="F")
AIC(mgcv.wus.null, mgcv.wus.tpa)

anova(mgcv.wus.t, mgcv.wus.tpa,test="F")
AIC(mgcv.wus.t, mgcv.wus.tpa)

# interaction
mgcv.wus.tp<- gam(influx ~ te(tann, pdsi, bs="tp", fx=T, k=4), data=wus )
summary(mgcv.wus.tp)
par(mfcol=c(1,1))
plot.gam(mgcv.wus.tp, residuals=T, se=T, shade=T, pch=16, cex=0.5)
vis.gam(mgcv.wus, plot.type="persp", type="response", color="gray", theta=-60)
gam.check(mgcv.wus.tp)
wus.fit <- predict(mgcv.wus.tp, wus)
wusin.fit <- predict(mgcv.wus.tp, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, mgcv.wus.tp$fitted.values); r2 <- r^2; r; r2
anova(mgcv.wus.null, mgcv.wus.tp,test="F")
AIC(mgcv.wus.null, mgcv.wus.tp)

anova(mgcv.wus.t, mgcv.wus.tp,test="F")
AIC(mgcv.wus.t, mgcv.wus.tp)

# tann + interaction
mgcv.wus.ttp <- gam(influx ~ te(tann, bs="tp", fx=T, k=4) + te(tann, pdsi, bs="tp", fx=F, k=4), data=wus )
summary(mgcv.wus.ttp)
par(mfcol=c(1,1))
plot.gam(mgcv.wus.ttp, residuals=T, se=T, shade=T, pch=16, cex=0.5)
vis.gam(mgcv.wus, plot.type="persp", type="response", color="gray", theta=-60)
gam.check(mgcv.wus.ttp)
wus.fit <- predict(mgcv.wus.ttp, wus)
wusin.fit <- predict(mgcv.wus.ttp, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, mgcv.wus.ttp$fitted.values); r2 <- r^2; r; r2
anova(mgcv.wus.null, mgcv.wus.ttp,test="F")
AIC(mgcv.wus.null, mgcv.wus.ttp)

anova(mgcv.wus.tp, mgcv.wus.ttp,test="F")
AIC(mgcv.wus.tp, mgcv.wus.ttp)

# pdsi + interaction
mgcv.uws.ptp <- gam(influx ~ te(pdsi, bs="tp", k=2) + te(tann, pdsi, bs="tp", k=4), data=wus )
summary(mgcv.uws.ptp)
par(mfcol=c(1,1))
plot.gam(mgcv.uws.ptp, residuals=T, se=T, shade=T, pch=16, cex=0.5)
vis.gam(mgcv.wus, plot.type="persp", type="response", color="gray", theta=-60)
gam.check(mgcv.uws.ptp)
wus.fit <- predict(mgcv.uws.ptp, wus)
wusin.fit <- predict(mgcv.uws.ptp, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, mgcv.uws.ptp$fitted.values); r2 <- r^2; r; r2
anova(mgcv.wus.null, mgcv.uws.ptp,test="F")
AIC(mgcv.wus.null, mgcv.uws.ptp)

anova(mgcv.wus.tp, mgcv.uws.ptp,test="F")
AIC(mgcv.wus.tp, mgcv.uws.ptp)

anova(mgcv.wus.ttp, mgcv.uws.ptp,test="F")
AIC(mgcv.wus.ttp, mgcv.uws.ptp)

# full model
mgcv.wus <- gam(influx ~ te(tann, bs="tp", fx=T, k=2) + te(pdsi, bs="tp", fx=T, k=2) + te(tann, pdsi, bs="tp", fx=F, k=2), data=wus )

summary(mgcv.wus)
par(mfcol=c(1,1))
plot.gam(mgcv.wus, residuals=T, se=T, shade=T, pch=16, cex=0.5)
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

anova(mgcv.wus.ttp, mgcv.wus,test="F")
AIC(mgcv.wus.ttp, mgcv.wus)

# plots
new.tann <- seq(from=-0.3, to=0.15, len=40); new.pdsi <- seq(from=-0.4, to=0.0, len=40)
new.x <- expand.grid(tann=new.tann, pdsi=new.pdsi)
mgcv.sfc <- predict(mgcv.wus, new.x)
#mgcv.sfc[new.x$pdsi >= 0.39] <- NA
#mgcv.sfc[new.x$pdsi >= -0.20 & new.x$tann >= -0.05] <- NA
#mgcv.sfc[new.x$pdsi >= 0.365 & new.x$tann <= -0.1] <- NA
#mgcv.sfc[new.x$pdsi <= 0.31 & new.x$tann >= -0.075] <- NA
mgcv.sfc <- matrix(mgcv.sfc, 40, 40)

par(mfcol=c(1,1))
persp(new.tann, new.pdsi, mgcv.sfc, theta=-50, phi=25, d=1.5, col="gray",ticktype="detailed",
      zlim=c(-0.5,0.5),
      xlab="MAT (C)", ylab="PDSI (JJA)",
      zlab = "Z-Score of Transformed Charcoal Influx")

image(new.tann, new.pdsi, mgcv.sfc)
points(wus$tann, wus$pdsi)

wus.fit <- predict(mgcv.wus, wus)
wusin.fit <- predict(mgcv.wus, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")

out.wus.fit <- cbind(wusin$Year,wusin$influx,wusin.fit)
colnames(out.wus.fit) <- c("Year","Influx","Pred")
write.table(out.wus.fit, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")