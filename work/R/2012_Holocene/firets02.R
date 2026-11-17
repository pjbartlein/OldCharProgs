dataname <- "e:\\Projects\\firets\\data\\euro_influx_tann.csv"
fts <- read.csv(dataname)
attach(fts)

influx.ts <- as.ts(influx, start=agets[1], end=agets[length(agets)], deltat=20)
tann.ts <- as.ts(tann, start=agets[1], end=agets[length(agets)], deltat=20)


ar.mod01 <- ar(ts.union(influx.ts,tann.ts), order.max=2)
acf(ar.mod01$resid[,1], na.action=na.pass)
acf(ar.mod01$resid[,2], na.action=na.pass)

arima.mod01 <- arima(influx.ts, order=c(1,1,1))
arima.mod01
acf(arima.mod01$resid, na.action=na.pass)

arima.mod02 <- arima(influx.ts, order=c(1,1,1), xreg=tann.ts)
arima.mod02
acf(arima.mod02$resid, na.action=na.pass)


lm.01 <- lm(influx.ts ~ tann.ts)
summary(lm.01)
AIC(lm.01)
acf(lm.01$resid)
Box.test(lm.01$resid, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 1)

library(nlme)

gls.01 <- gls(influx.ts ~ tann.ts, correlation=corARMA(p=0,q=1))
summary(gls.01)
AIC(gls.01)
acf(gls.01$resid)
Box.test(gls.01$resid, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 2)

arima.tann.01 <- arima(tann.ts, order=c(1,0,1))
arima.tann.01
AIC(arima.tann.01)
acf(arima.tann.01$resid)
Box.test(arima.tann.01$resid, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 2)
arima.tann.01$fit <- tann.ts + arima.tann.01$resid
plot(tann.ts ~ agets, pch=16, cex=0.5)
lines(arima.tann.01$fit ~ agets, col="red")
plot(tann.ts ~ arima.tann.01$fit)
plot(arima.tann.01$resid)
