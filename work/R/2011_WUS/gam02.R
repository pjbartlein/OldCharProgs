
library(mgcv)

#dataname <- "e:\\Projects\\GPWG\\work\\QSR\\firets\\data\\euro_influx_tann_200.csv"
dataname <- "e:\\Projects\\GPWG\\work\\QSR\\firets\\data\\namer_influx_tann_200.csv"
fts <- read.csv(dataname)
head(fts)

endage <- max(fts$age)
begage <- min(fts$age)
delta <- fts$age[2]-fts$age[1]

# convert data to time series class
age.ts <- ts(fts$age, start=begage,  end=endage, deltat=delta)
influx.ts <- ts(fts$influx, start=begage,  end=endage, deltat=delta)
tann.ts <- ts(fts$tann, start=begage,  end=endage, deltat=delta)

# plot the data, and get ACF and PACF and Ljung-Box statistic
plot(influx.ts, type="o")
acf(influx.ts)
pacf(influx.ts)
Box.test(influx.ts, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 3)

# fit a GAM to influx data using tann as a predictor 
influx.gam01 <- gam(influx.ts ~ te(tann.ts, k=4))
summary(influx.gam01)
AIC(influx.gam01)
plot(influx.gam01)

influx.gam01.fit <- ts(influx.gam01$fitted.values, start=begage,  end=endage, deltat=delta)
influx.gam01.resid <- ts(influx.gam01$residuals, start=begage,  end=endage, deltat=delta)
plot(influx.ts, type="o")
lines(influx.gam01.fit, col="orange")

# plot the residuals and get ACF, PACF and Ljung-Box statistic
ts.plot(influx.gam01.resid)
acf(influx.gam01.resid)
pacf(influx.gam01.resid)
Box.test(influx.gam01.resid, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 2)

# now ARMAX -- skip

influx.arima01 <- arima(influx.ts, xreg=influx.gam01.fit, order=c(0,1,1)) # tann.ts  tann.lm01.fit
influx.arima01
AIC(influx.arima01)
influx.arima01.resid <- ts(influx.arima01$resid, start=begage,  end=endage, deltat=delta)
influx.arima01.fit <- ts(influx.ts-influx.arima01$resid, start=begage,  end=endage, deltat=delta)
plot(influx.ts, type="o")
lines(influx.arima01.fit, col="green")

plot(influx.arima01.resid)
acf(influx.arima01.resid)
Box.test(influx.arima01$resid, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 3)

plot(influx.ts, type="o")
lines(influx.arima01.fit, col="purple")

R <- cor(influx.ts, influx.arima01.fit)
Rsq <- R^2
print (c(R,Rsq))

# ARIMA of GAM residuals

resid.arima01 <- arima(influx.gam01.resid,  order=c(0,1,1)) 
resid.arima01
AIC(resid.arima01)
resid.arima01.resid <- ts(resid.arima01$resid, start=begage,  end=endage, deltat=delta)
resid.arima01.fit <- ts(influx.gam01.resid-resid.arima01$resid, start=begage,  end=endage, deltat=delta)
plot(influx.gam01.resid, type="o")
lines(resid.arima01.fit, col="green")

plot(resid.arima01.resid)
acf(resid.arima01.resid)
Box.test(resid.arima01.resid, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 3)

plot(influx.ts, type="o")
lines(influx.gam01.fit+resid.arima01.fit, col="purple")

R <- cor(influx.ts, influx.gam01.fit+resid.arima01.fit)
Rsq <- R^2
print (c(R,Rsq))


# Europe ARIMA(1,0,1)
# North America ARIMA(0,1,1)
