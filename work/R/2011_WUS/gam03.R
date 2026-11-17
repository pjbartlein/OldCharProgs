# GAM plus ARIMA models of influx, North America and Europe
library(mgcv)

datafile <- "e:\\Projects\\GPWG\\work\\QSR\\firets\\data\\euro_influx_tann_200.csv"
outfile <- "e:\\Projects\\GPWG\\work\\QSR\\firets\\data\\euro_200out.csv"
#datafile <- "e:\\Projects\\GPWG\\work\\QSR\\firets\\data\\namer_influx_tann_200.csv"
#outfile <- "e:\\Projects\\GPWG\\work\\QSR\\firets\\data\\namer_200_out.csv"
fts <- read.csv(datafile)
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
plot(influx.gam01, residuals=TRUE)

# plot the observed and fitted influx values
influx.gam01.fit <- ts(influx.gam01$fitted.values, start=begage,  end=endage, deltat=delta)
influx.gam01.resid <- ts(influx.gam01$residuals, start=begage,  end=endage, deltat=delta)
plot(influx.ts, type="o")
lines(influx.gam01.fit, col="orange", lwd=2)

# plot the residuals, and get ACF, PACF and Ljung-Box statistic
ts.plot(influx.gam01.resid)
acf(influx.gam01.resid)
pacf(influx.gam01.resid)
Box.test(influx.gam01.resid, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 2)


# ARIMA of GAM residuals (influx.gam01.resid)

# ARMA(1,1) first (final for Europe)
resid.arima01 <- arima(influx.gam01.resid,  order=c(1,0,1)) 
resid.arima01
AIC(resid.arima01)

# ARIMA(1,0,1) (final for North America)
resid.arima01 <- arima(influx.gam01.resid,  order=c(0,1,1)) 
resid.arima01
AIC(resid.arima01)

# plot the observed and fitted GAM residual series
resid.arima01.resid <- ts(resid.arima01$resid, start=begage,  end=endage, deltat=delta)
resid.arima01.fit <- ts(influx.gam01.resid-resid.arima01$resid, start=begage,  end=endage, deltat=delta)
plot(influx.gam01.resid, type="o")
lines(resid.arima01.fit, col="green", lwd=2)

# plot the residuals from the ARIMA model, and get ACF, PACF and Ljung-Box statistic
plot(resid.arima01.resid)
acf(resid.arima01.resid)
Box.test(resid.arima01.resid, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 3)

# plot the influx data and combined fitted values
influx.fit <- influx.gam01.fit+resid.arima01.fit
plot(influx.ts, type="o")
lines(influx.fit, col="purple", lwd=2)

# poor man's R^2
R <- cor(influx.ts, influx.gam01.fit+resid.arima01.fit)
Rsq <- R^2
print (c(R,Rsq))

# write out some data
out.fit <- cbind(fts$agets,tann.ts,influx.ts,influx.gam01.fit,influx.gam01.resid,
                 resid.arima01.fit,resid.arima01.resid,influx.fit) 
#colnames(out.wus.fit) <- c("Age","tann","influx")
write.table(out.fit, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")


# Europe ARIMA(1,0,1)
# North America ARIMA(0,1,1)

