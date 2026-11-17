# wus influx
# library(mgcv)

workpath <-"e:\\Projects\\GPWG\\work\\wus\\R\\"
curvesname <- "boot100.csv"
outname <- "boot100_fitted_1800.csv"
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
mgcv.wus.t <- gam(influx ~ te(tann, bs="tp", k=4), data=wus )
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
mgcv.wus.d <- gam(influx ~ te(dai, bs="tp", k=4), data=wus )
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

# interaction
mgcv.wus.td <- gam(influx ~ te(tann, dai, bs="tp", k=4), data=wus )
summary(mgcv.wus.td)
plot.gam(mgcv.wus.td, residuals=T, se=T, shade=T, pch=16, cex=0.5)
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
mgcv.wus.ttd <- gam(influx ~ te(tann, bs="tp", k=2) + te(tann, dai, bs="tp", k=4), data=wus )
summary(mgcv.wus.ttd)
plot.gam(mgcv.wus.ttd, residuals=T, se=T, shade=T, pch=16, cex=0.5)
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
lm.uws.dtd <- gam(influx ~ te(dai, bs="tp", k=2) + te(tann, dai, bs="tp", k=4), data=wus )
summary(lm.uws.dtd)
plot.gam(lm.uws.dtd, residuals=T, se=T, shade=T, pch=16, cex=0.5)
gam.check(lm.uws.dtd)
wus.fit <- predict(lm.uws.dtd, wus)
wusin.fit <- predict(lm.uws.dtd, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, lm.uws.dtd$fitted.values); r2 <- r^2; r; r2
anova(mgcv.wus.null, lm.uws.dtd,test="F")
AIC(mgcv.wus.null, lm.uws.dtd)

anova(mgcv.wus.td, lm.uws.dtd,test="F")
AIC(mgcv.wus.td, lm.uws.dtd)

anova(mgcv.wus.ttd, lm.uws.dtd,test="F")
AIC(mgcv.wus.ttd, lm.uws.dtd)

# full model
mgcv.wus <- gam(influx ~ te(tann, bs="tp", k=2) + te(dai, bs="tp", k=2) + te(tann, dai, bs="tp", k=2), data=wus )
summary(mgcv.wus)
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

anova(mgcv.wus.ttd, mgcv.wus,test="F")
AIC(mgcv.wus.ttd, mgcv.wus)

# plots
new.tann <- seq(from=-0.3, to=0.1, len=40); new.pdsi <- seq(from=-0.3, to=0.0, len=38)
new.dai <- seq(from=0.3, to=0.40, len=40)
##new.wus$tmp <- seq(from=min(wus$tmp), to=max(wus$tmp), len=40); new.wus$pme <- seq(from=min(wus$pme), to=max(wus$pme), len=40)
new.x <- expand.grid(tann=new.tann, dai=new.dai)
lm.sfc <- predict(mgcv.wus, new.x)
#lm.sfc[new.x$dai >= 0.39] <- NA
#lm.sfc[new.x$dai >= 0.375 & new.x$tann <= -0.05] <- NA
#lm.sfc[new.x$dai >= 0.365 & new.x$tann <= -0.1] <- NA
#lm.sfc[new.x$dai <= 0.31 & new.x$tann >= -0.075] <- NA
lm.sfc <- matrix(lm.sfc, 40, 40)

par(mfcol=c(1,1))
persp(new.tann, new.dai, lm.sfc, theta=-50, phi=30, d=1.5, col="gray",ticktype="detailed",
      zlim=c(-0.5,0.5),
      xlab="MAT (C)", ylab="DAI (JJA)",
      zlab = "Z-Score of Transformed Charcoal Influx")

image(new.tann, new.dai, lm.sfc)
points(wus$tann, wus$dai)

wus.fit <- predict(mgcv.wus, wus)
wusin.fit <- predict(mgcv.wus, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")

out.wus.fit <- cbind(wusin$Year,wusin$influx,wusin.fit)
colnames(out.wus.fit) <- c("Year","Influx","Pred")
write.table(out.wus.fit, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")