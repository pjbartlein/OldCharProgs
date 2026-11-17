library(nlme)

dataname <- "e:\\Projects\\GPWG\\work\\QSR\\firets\\data\\euro_influx_tann_200.csv"
fts <- read.csv(dataname)
#attach(fts)
head(fts)

endage <- max(fts$age)
begage <- min(fts$age)
delta <- fts$age[2]-fts$age[1]

age.ts <- ts(fts$age, start=begage,  end=endage, deltat=delta)
influx.ts <- ts(fts$influx, start=begage,  end=endage, deltat=delta)
tann.ts <- ts(fts$tann, start=begage,  end=endage, deltat=delta)

tann.lm01 <- lm(tann.ts ~ age.ts + I(age.ts^2) + I(age.ts^3))
summary(tann.lm01)
AIC(tann.lm01)
tann.lm01.fit <- ts(tann.lm01$fitted.values, start=begage,  end=endage, deltat=delta)
tann.lm01.resid <- ts(tann.lm01$residuals, start=begage,  end=endage, deltat=delta)
ts.plot(tann.ts)
lines(tann.lm01.fit, col="red")
ts.plot(tann.lm01.resid)

influx.lm01 <- lm(influx.ts ~ tann.lm01.fit + tann.lm01.resid)
summary(influx.lm01)
AIC(influx.lm01)
influx.lm01.fit <- ts(influx.lm01$fitted.values, start=begage,  end=endage, deltat=delta)
influx.lm01.resid <- ts(influx.lm01$residuals, start=begage,  end=endage, deltat=delta)
ts.plot(influx.ts)
lines(influx.lm01.fit, col="red")

ts.plot(influx.lm01.resid)
acf(influx.lm01.resid)
pacf(influx.lm01.resid)

glsControl(msMaxIter = 20, msVerbose = TRUE, singular.ok = TRUE)
influx.gls01 <- gls(influx.ts ~ tann.lm01.fit + tann.lm01.resid, correlation=corARMA(p=2, q=0))
summary(influx.gls01)
AIC(influx.gls01)
acf(influx.gls01$resid)
Box.test(influx.gls01$resid, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 2)




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

ic <- read.csv(file="e:\\Projects\\paleovar\\data\\datasets\\ic2_20yr.csv")  

endage <- max(ic$Age)
begage <- min(ic$Age)
delta <- ic$Age[2]-ic$Age[1]
age=ts(ic$AgeBP1950, start=begage, end=endage, deltat=delta)
d18O=ts(ic$NGRIP_d18O, start=begage, end=endage, deltat=delta)
ngrip <- ts.intersect(age, d18O)

time(age)