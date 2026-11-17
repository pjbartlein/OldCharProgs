# GAM plus ARIMA models of influx, North America
library(mgcv)

datafile <- "e:\\Projects\\GPWG\\work\\QSR\\firets\\data\\euro_influx_tann_200.csv"
outfile <- "e:\\Projects\\GPWG\\work\\QSR\\firets\\data\\euro_200_out.csv"
gamplotfile <- "e:\\Projects\\GPWG\\work\\QSR\\firets\\data\\euro_200_gamplot.csv"
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
Box.test(influx.ts.euro, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 3)

# fit a GAM to influx data using tann as a predictor 
influx.gam01.euro <- gam(influx.ts.euro ~ te(tann.ts.euro, k=4))
summary(influx.gam01.euro)
AIC(influx.gam01.euro)
plot(influx.gam01.euro, residuals=TRUE, pch=16, cex=1.0)

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
#AIC(resid.arima01.euro)

# plot the observed and fitted GAM residual series
resid.arima01.euro.resid <- ts(resid.arima01.euro$resid, start=begage,  end=endage, deltat=delta)
resid.arima01.euro.fit <- ts(influx.gam01.euro.resid-resid.arima01.euro$resid, start=begage,  end=endage, deltat=delta)
plot(influx.gam01.euro.resid, type="o")
lines(resid.arima01.euro.fit, col="green", lwd=2)

# plot the residuals from the ARIMA model, and get ACF, PACF and Ljung-Box statistic
plot(resid.arima01.euro.resid)
acf(resid.arima01.euro.resid)
Box.test(resid.arima01.euro.resid, lag = 10, type = c("Box-Pierce", "Ljung-Box"), fitdf = 3)

# plot the influx data and combined fitted values
influx.fit.euro <- influx.gam01.euro.fit+resid.arima01.euro.fit
plot(influx.ts.euro, type="o")
lines(influx.fit.euro, col="purple", lwd=2)

# poor man's R^2
R <- cor(influx.ts.euro, influx.gam01.euro.fit+resid.arima01.euro.fit)
Rsq <- R^2
print (c(R,Rsq))

# write out some data
out.fit <- cbind(fts.euro$agets,tann.ts.euro,influx.ts.euro,influx.gam01.euro.fit,influx.gam01.euro.resid,
                 resid.arima01.euro.fit,resid.arima01.euro.resid,influx.fit.euro) 
#colnames(out.wus.fit) <- c("Age","tann","influx")
write.table(out.fit, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")

# extract components of GAM plot
fv <- predict(influx.gam01.euro,type="terms", se.fit=TRUE) ## get term estimates       
partial.resid.euro <- residuals(influx.gam01.euro,type="working") + fv$fit 
se.fit.euro <- fv$se.fit
se.upper.euro <- fv$fit + 2.0*fv$se.fit
se.lower.euro <- fv$fit - 2.0*fv$se.fit

# check extraction by overplotting
plot(influx.gam01.euro, ylim=c(-0.8,0.4), residuals=TRUE, pch=16 ) # GAM plot
lines(fts.euro$tann, fv$fit, col="red", lwd=0.2)
lines(fts.euro$tann, se.upper.euro, col="red", lwd=0.2)
lines(fts.euro$tann, se.lower.euro, col="red", lwd=0.2)
points(tann.ts.euro, partial.resid.euro, pch=16, col="red", cex=0.4) # partial residuals

# write GAM plot components
gamplt <- cbind(fts.euro$tann,fv$fit,se.upper.euro,se.lower.euro,partial.resid.euro)          
colnames(gamplt) <- c("euro.tann","euro.tann_fit","euro.tann.upper",
                      "euro.tann.lower","euro.tann.partialresid")
write.table(gamplt, file=gamplotfile, row.names=FALSE, col.names=TRUE, sep=",")