# wus influx
# library(mgcv)

workpath <-"e:\\Projects\\GPWG\\work\\wus\\R\\"
curvesname <- "boot100.csv"
outname <- "boot100_fitted.csv"
outfile <- paste(workpath, outname, sep="")
                 
wusin <- read.csv(paste(workpath, curvesname, sep=""))
summary(wusin)
wus <- wusin[wusin$Year <= 1850,]
#wus <- wusin
summary(wus)

par(mfcol=c(1,1))

# regressions autocorrelated residuals

mgcv.null.wus <- gam(influx ~ 1, data=wus)
summary(mgcv.null.wus)

#mgcv.wus <- gam(influx ~ te(tann, bs="tp", k=2) + te(dai, bs="tp", k=2) + te(tann, dai, bs="tp", k=4) , data=wus )

# tann ARMA residuals
mgcv.wus.t <- gamm(influx ~ te(tann, bs="tp", k=2) + te(dai, bs="tp", k=4) + te(tann, dai, bs="tp", k=2),
                  correlation=corARMA(p=0, q=1), data=wus )
summary(mgcv.wus.t$gam); summary(mgcv.wus.t$lme)

plot.gam(mgcv.wus.t$gam, residuals=T, se=T, shade=T, pch=16, cex=0.5)
gam.check(mgcv.wus.t$gam)
vis.gam(mgcv.wus.t$gam)
acf(mgcv.wus.t$lme$residuals)

wus.fit <- predict(mgcv.wus.t$gam, wus)
wusin.fit <- predict(mgcv.wus.t$gam, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, mgcv.wus.t$gam$fitted.values); r2 <- r^2; r; r2



