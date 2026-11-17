
library(mgcv)

dataname <- "e:\\Projects\\GPWG\\work\\QSR\\firets\\data\\euro_influx_tann_200.csv"
#dataname <- "e:\\Projects\\GPWG\\work\\QSR\\firets\\data\\namer_influx_tann_200.csv"
fts <- read.csv(dataname)
attach(fts)
head(fts)

endage <- max(fts$age)
begage <- min(fts$age)
delta <- fts$age[2]-fts$age[1]

age.ts <- ts(fts$age, start=begage,  end=endage, deltat=delta)
influx.ts <- ts(fts$influx, start=begage,  end=endage, deltat=delta)
tann.ts <- ts(fts$tann, start=begage,  end=endage, deltat=delta)

acf(influx.ts)
pacf(influx.ts)
Box.test(influx.ts, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 3)

# GAM 
influx.gam01 <- gam(influx.ts ~ s(tann.ts)) #, correlation=corARMA(p=1, q=1))
summary(influx.gam01)
AIC(influx.gam01)
plot(influx.gam01)

influx.gam01.fit <- ts(influx.gam01$fitted.values, start=begage,  end=endage, deltat=delta)
influx.gam01.resid <- ts(influx.gam01$residuals, start=begage,  end=endage, deltat=delta)
plot(influx.ts, type="o")
lines(influx.gam01.fit, col="orange")

ts.plot(influx.gam01.resid)
acf(influx.gam01.resid)
pacf(influx.gam01.resid)
Box.test(influx.gam01.resid, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 2)

# GAM with autocorrelated errors
influx.gam02<- gamm(influx ~ te(tann, k=10) + te(agets, k=10), correlation=corARMA(p=1, q=0))
summary(influx.gam02$gam)
intervals(influx.gam02$lme)
#AIC(influx.gam02$gam)
plot(influx.gam02$gam)

influx.gam02.fit <- ts(influx.gam02$gam$fitted.values, start=begage,  end=endage, deltat=delta)
influx.gam02.resid <- ts(influx.gam02$gam$residuals, start=begage,  end=endage, deltat=delta)
plot(influx.ts, type="o")
lines(influx.gam02.fit, col="orange")

ts.plot(influx.gam01.resid)
acf(influx.gam02.resid)
pacf(influx.gam02.resid)
Box.test(influx.gam02.resid, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 2)