# GAM plus ARIMA models of influx, North America and Europe
library(mgcv)

datafile <- "e:\\Projects\\GPWG\\work\\QSR\\firets\\data\\namer_influx_tann_200.csv"
outfile <- "e:\\Projects\\GPWG\\work\\QSR\\firets\\data\\namer_200_out.csv"
fts.namer <- read.csv(datafile)
head(fts.namer)

endage <- max(fts.namer$age)
begage <- min(fts.namer$age)
delta <- fts.namer$age[2]-fts.namer$age[1]

# convert data to time series class
age.ts.namer <- ts(fts.namer$age, start=begage,  end=endage, deltat=delta)
influx.ts.namer <- ts(fts.namer$influx, start=begage,  end=endage, deltat=delta)
tann.ts.namer <- ts(fts.namer$tann, start=begage,  end=endage, deltat=delta)

# plot the data, and get ACF and PACF and Ljung-Box statistic
plot(influx.ts.namer, type="o")
acf(influx.ts.namer)
pacf(influx.ts.namer)
Box.test(influx.ts, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 3)

# fit a GAM to influx data using tann as a predictor 
influx.gam01.namer <- gam(influx.ts.namer ~ te(tann.ts.namer, k=4))
summary(influx.gam01.namer)
AIC(influx.gam01.namer)
plot(influx.gam01.namer, residuals=TRUE, pch=16, cex=1.0)

# plot the observed and fitted influx values
influx.gam01.namer.fit <- ts(influx.gam01.namer$fitted.values, start=begage,  end=endage, deltat=delta)
influx.gam01.namer.resid <- ts(influx.gam01.namer$residuals, start=begage,  end=endage, deltat=delta)
plot(influx.ts.namer, type="o")
lines(influx.gam01.namer.fit, col="orange", lwd=2)

# plot the residuals, and get ACF, PACF and Ljung-Box statistic
ts.plot(influx.gam01.namer.resid)
acf(influx.gam01.namer.resid)
pacf(influx.gam01.namer.resid)
Box.test(influx.gam01.namer.resid, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 2)


# ARIMA of GAM residuals (influx.gam01.resid)

# ARMA(1,1) first (final for Europe)
#resid.arima01 <- arima(influx.gam01.resid,  order=c(1,0,1)) 
#resid.arima01
#AIC(resid.arima01)

# ARIMA(1,0,1) (final for North America)
resid.arima01.namer <- arima(influx.gam01.namer.resid,  order=c(0,1,1)) 
resid.arima01.namer
AIC(resid.arima01)

# plot the observed and fitted GAM residual series
resid.arima01.namer.resid <- ts(resid.arima01.namer$resid, start=begage,  end=endage, deltat=delta)
resid.arima01.namer.fit <- ts(influx.gam01.namer.resid-resid.arima01.namer$resid, start=begage,  end=endage, deltat=delta)
plot(influx.gam01.namer.resid, type="o")
lines(resid.arima01.namer.fit, col="green", lwd=2)

# plot the residuals from the ARIMA model, and get ACF, PACF and Ljung-Box statistic
plot(resid.arima01.namer.resid)
acf(resid.arima01.namer.resid)
Box.test(resid.arima01.namer.resid, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 3)

# plot the influx data and combined fitted values
influx.fit.namer <- influx.gam01.namer.fit+resid.arima01.namer.fit
plot(influx.ts.namer, type="o")
lines(influx.fit.namer, col="purple", lwd=2)

# poor man's R^2
R <- cor(influx.ts.namer, influx.gam01.namer.fit+resid.arima01.namer.fit)
Rsq <- R^2
print (c(R,Rsq))

# write out some data
out.fit <- cbind(fts$agets.namer,tann.ts.namer,influx.ts.namer,influx.gam01.namer.fit,influx.gam01.namer.resid,
                 resid.arima01.fit,resid.arima01.resid,influx.fit) 
#colnames(out.wus.fit) <- c("Age","tann","influx")
write.table(out.fit, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")


datafile <- "e:\\Projects\\GPWG\\work\\QSR\\firets\\data\\euro_influx_tann_200.csv"
outfile <- "e:\\Projects\\GPWG\\work\\QSR\\firets\\data\\euro_200_out.csv"

fts.euro <- read.csv(datafile)
head(fts.euro)

endage <- max(fts.euro$age)
begage <- min(fts.euro$age)
delta <- fts.euro$age[2]-fts.euro$age[1]

# convert data to time series class
age.ts.euro <- ts(fts.euro$age, start=begage,  end=endage, deltat=delta)
influx.ts.euro <- ts(fts.euro$influx, start=begage,  end=endage, deltat=delta)
tann.ts.euro <- ts(fts.euro$tann, start=begage,  end=endage, deltat=delta)

# plot the data, and get ACF and PACF and Ljung-Box statistic
plot(influx.ts.euro, type="o")
acf(influx.ts.euro)
pacf(influx.ts.euro)
Box.test(influx.ts, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 3)

# fit a GAM to influx data using tann as a predictor 
influx.gam01.euro <- gam(influx.ts.euro ~ te(tann.ts.euro, k=4))
summary(influx.gam01.euro)
AIC(influx.gam01.euro)
plot(influx.gam01.euro, residuals=TRUE, pch=16)

# plot the observed and fitted influx values
influx.gam01.euro.fit <- ts(influx.gam01.euro$fitted.values, start=begage,  end=endage, deltat=delta)
influx.gam01.euro.resid <- ts(influx.gam01.euro$residuals, start=begage,  end=endage, deltat=delta)
plot(influx.ts.euro, type="o")
lines(influx.gam01.euro.fit, col="orange", lwd=2)

# plot the residuals, and get ACF, PACF and Ljung-Box statistic
ts.plot(influx.gam01.euro.resid)
acf(influx.gam01.euro.resid)
pacf(influx.gam01.euro.resid)
Box.test(influx.gam01.euro.resid, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 2)


# ARIMA of GAM residuals (influx.gam01.resid)

# ARMA(1,1) first (final for Europe)
resid.arima01.euro <- arima(influx.gam01.euro.resid,  order=c(1,0,1)) 
resid.arima01.euro
AIC(resid.arima01.euro)

# ARIMA(1,0,1) (final for North America)
#resid.arima01.euro <- arima(influx.gam01.euro.resid,  order=c(0,1,1)) 
#resid.arima01.euro
#AIC(resid.arima01)

# plot the observed and fitted GAM residual series
resid.arima01.euro.resid <- ts(resid.arima01.euro$resid, start=begage,  end=endage, deltat=delta)
resid.arima01.euro.fit <- ts(influx.gam01.euro.resid-resid.arima01.euro$resid, start=begage,  end=endage, deltat=delta)
plot(influx.gam01.euro.resid, type="o")
lines(resid.arima01.euro.fit, col="green", lwd=2)

# plot the residuals from the ARIMA model, and get ACF, PACF and Ljung-Box statistic
plot(resid.arima01.euro.resid)
acf(resid.arima01.euro.resid)
Box.test(resid.arima01.euro.resid, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 3)

# plot the influx data and combined fitted values
influx.fit.euro <- influx.gam01.euro.fit+resid.arima01.euro.fit
plot(influx.ts.euro, type="o")
lines(influx.fit.euro, col="purple", lwd=2)

# poor man's R^2
R <- cor(influx.ts.euro, influx.gam01.euro.fit+resid.arima01.euro.fit)
Rsq <- R^2
print (c(R,Rsq))

# write out some data
out.fit <- cbind(fts$agets.euro,tann.ts.euro,influx.ts.euro,influx.gam01.euro.fit,influx.gam01.euro.resid,
                 resid.arima01.fit,resid.arima01.resid,influx.fit) 
#colnames(out.wus.fit) <- c("Age","tann","influx")
write.table(out.fit, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")

datafile <- "e:\\Projects\\GPWG\\work\\QSR\\firets\\data\\asia_influx_d18O_200.csv"
outfile <- "e:\\Projects\\GPWG\\work\\QSR\\firets\\data\\asia_200_out.csv"

fts.asia <- read.csv(datafile)
head(fts.asia)

endage <- max(fts.asia$age)
begage <- min(fts.asia$age)
delta <- fts.asia$age[2]-fts.asia$age[1]

# convert data to time series class
age.ts.asia <- ts(fts.asia$age, start=begage,  end=endage, deltat=delta)
influx.ts.asia <- ts(fts.asia$influx, start=begage,  end=endage, deltat=delta)
d18O.ts.asia <- ts(fts.asia$d18O, start=begage,  end=endage, deltat=delta)

# plot the data, and get ACF and PACF and Ljung-Box statistic
plot(influx.ts.asia, type="o")
acf(influx.ts.asia)
pacf(influx.ts.asia)
Box.test(influx.ts, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 3)

# fit a GAM to influx data using d18O as a predictor 
influx.gam01.asia <- gam(influx.ts.asia ~ te(d18O.ts.asia, bs="tp", k=5))
summary(influx.gam01.asia)
AIC(influx.gam01.asia)
plot(influx.gam01.asia, residuals=TRUE, pch=16)

# plot the observed and fitted influx values
influx.gam01.asia.fit <- ts(influx.gam01.asia$fitted.values, start=begage,  end=endage, deltat=delta)
influx.gam01.asia.resid <- ts(influx.gam01.asia$residuals, start=begage,  end=endage, deltat=delta)
plot(influx.ts.asia, type="o")
lines(influx.gam01.asia.fit, col="orange", lwd=2)

# plot the residuals, and get ACF, PACF and Ljung-Box statistic
ts.plot(influx.gam01.asia.resid)
acf(influx.gam01.asia.resid)
pacf(influx.gam01.asia.resid)
Box.test(influx.gam01.asia.resid, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 2)


# ARIMA of GAM residuals (influx.gam01.resid)

# ARMA(1,1) first (final for Europe)
resid.arima01.asia <- arima(influx.gam01.asia.resid,  order=c(1,0,1)) 
resid.arima01.asia
AIC(resid.arima01.asia)

# ARIMA(1,0,1) (final for North America)
#resid.arima01.asia <- arima(influx.gam01.asia.resid,  order=c(0,1,1)) 
#resid.arima01.asia
#AIC(resid.arima01)

# plot the observed and fitted GAM residual series
resid.arima01.asia.resid <- ts(resid.arima01.asia$resid, start=begage,  end=endage, deltat=delta)
resid.arima01.asia.fit <- ts(influx.gam01.asia.resid-resid.arima01.asia$resid, start=begage,  end=endage, deltat=delta)
plot(influx.gam01.asia.resid, type="o")
lines(resid.arima01.asia.fit, col="green", lwd=2)

# plot the residuals from the ARIMA model, and get ACF, PACF and Ljung-Box statistic
plot(resid.arima01.asia.resid)
acf(resid.arima01.asia.resid)
Box.test(resid.arima01.asia.resid, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 3)

# plot the influx data and combined fitted values
influx.fit.asia <- influx.gam01.asia.fit+resid.arima01.asia.fit
plot(influx.ts.asia, type="o")
lines(influx.fit.asia, col="purple", lwd=2)

# poor man's R^2
R <- cor(influx.ts.asia, influx.gam01.asia.fit+resid.arima01.asia.fit)
Rsq <- R^2
print (c(R,Rsq))

# write out some data
out.fit <- cbind(fts$agets.asia,d18O.ts.asia,influx.ts.asia,influx.gam01.asia.fit,influx.gam01.asia.resid,
                 resid.arima01.fit,resid.arima01.resid,influx.fit) 
#colnames(out.wus.fit) <- c("Age","d18O","influx")
write.table(out.fit, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")


