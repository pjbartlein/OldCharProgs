
dataname <- "e:\\Projects\\GPWG\\work\\QSR\\firets\\data\\euro_influx_tann_200.csv"
#dataname <- "e:\\Projects\\GPWG\\work\\QSR\\firets\\data\\namer_influx_tann_200.csv"
fts <- read.csv(dataname)
#attach(fts)
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

tann.lm01 <- lm(tann.ts ~ age.ts + I(age.ts^2) + I(age.ts^3))
summary(tann.lm01)
AIC(tann.lm01)
tann.lm01.fit <- ts(tann.lm01$fitted.values, start=begage,  end=endage, deltat=delta)
tann.lm01.resid <- ts(tann.lm01$residuals, start=begage,  end=endage, deltat=delta)
plot(tann.ts, type="o")
lines(tann.lm01.fit, col="orange")
ts.plot(tann.lm01.resid)
acf(tann.lm01.resid)
Box.test(tann.lm01$resid, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 2)

#influx.lm01 <- lm(influx.ts ~ tann.lm01.fit + tann.lm01.resid )  # + tann.lm01.resid
influx.lm01 <- lm(influx.ts ~ tann.lm01.fit  + tann.lm01.resid  )  # 
summary(influx.lm01)
AIC(influx.lm01)
influx.lm01.fit <- ts(influx.lm01$fitted.values, start=begage,  end=endage, deltat=delta)
influx.lm01.resid <- ts(influx.lm01$residuals, start=begage,  end=endage, deltat=delta)
plot(influx.ts, type="o")
lines(influx.lm01.fit, col="red")
plot(influx.lm01.resid)
acf(influx.lm01.resid)
Box.test(influx.lm01$resid, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 2)

xm <- as.matrix(cbind(tann.lm01.fit,tann.lm01.resid))
influx.arima01 <- arima(influx.ts, xreg=xm, order=c(1,0,1)) # tann.ts  tann.lm01.fit
influx.arima01
AIC(influx.arima01)
influx.arima01.resid <- ts(influx.arima01$resid, start=begage,  end=endage, deltat=delta)
influx.arima01.fit <- ts(influx.ts-influx.arima01$resid, start=begage,  end=endage, deltat=delta)
plot(influx.ts, type="o")
lines(influx.arima01.fit, col="green")
plot(influx.arima01.resid)
acf(influx.arima01.resid)
Box.test(influx.arima01$resid, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 3)

#influx.sm.fit <- 0.1392 + 0.2862*tann.lm01.fit -0.0721*tann.lm01.fit^2
influx.sm.fit <- 0.1784 + 0.4856*tann.lm01.fit -0.0779*tann.lm01.resid
plot(influx.ts, type="o")
lines(influx.sm.fit, col="blue", lwd=2)
lines(influx.arima01.fit-influx.sm.fit, col="pink", lwd=2)

cor(influx.ts, influx.arima01.fit)

# best model for Europe and North America
#influx.lm01 <- lm(influx.ts ~ tann.lm01.fit + tann.lm01.resid )  
#influx.arima01 <- arima(influx.ts, xreg=tann.lm01.fit, order=c(1,0,1)) 
