# wus influx
# library(lm)

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

# regressions -- polynomials

lm.wus.null <- lm(influx ~ 1, data=wus)
summary(lm.wus.null)

#lm.wus <- gam(influx ~ te(tann, bs="tp", k=2) + te(dai, bs="tp", k=2) + te(tann, dai, bs="tp", k=4) , data=wus )

# tann alone
lm.wus.t <- lm(influx ~ tann +I(tann^2) +I(tann^2) + dai + I(dai^2) +I(tann*dai), data=wus )
summary(lm.wus.t)
par(mfcol=c(2,2))
plot(lm.wus.t, residuals=T, se=T, shade=T, pch=16, cex=0.5)

wus.fit <- predict(lm.wus.t, wus)
wusin.fit <- predict(lm.wus.t, wusin)
par(mfcol=c(1,1))
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, lm.wus.t$fitted.values); r2 <- r^2; r; r2
anova(lm.wus.null, lm.wus.t,test="F")
AIC(lm.wus.null, lm.wus.t)

