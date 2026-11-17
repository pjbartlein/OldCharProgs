# GAM plus ARIMA models of influx, North America
library(mgcv)

datafile <- "e:\\Projects\\GPWG\\work\\QSR\\firets\\data\\namer_influx_tann_200.csv"
outfile <- "e:\\Projects\\GPWG\\work\\QSR\\firets\\data\\namer_200_out.csv"
gamplotfile <- "e:\\Projects\\GPWG\\work\\QSR\\firets\\data\\namer_200_gamplot.csv"
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
Box.test(influx.ts.namer, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 3)

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
resid.arima01.namer <- arima(influx.gam01.namer.resid,  order=c(1,0,1)) 
resid.arima01.namer
AIC(resid.arima01.namer)

# ARIMA(1,0,1) (final for North America)
resid.arima01.namer <- arima(influx.gam01.namer.resid,  order=c(0,1,1)) 
resid.arima01.namer
AIC(resid.arima01.namer)

# plot the observed and fitted GAM residual series
resid.arima01.namer.resid <- ts(resid.arima01.namer$resid, start=begage,  end=endage, deltat=delta)
resid.arima01.namer.fit <- ts(influx.gam01.namer.resid-resid.arima01.namer$resid, start=begage,  end=endage, deltat=delta)
plot(influx.gam01.namer.resid, type="o")
lines(resid.arima01.namer.fit, col="green", lwd=2)

# plot the residuals from the ARIMA model, and get ACF, PACF and Ljung-Box statistic
plot(resid.arima01.namer.resid)
acf(resid.arima01.namer.resid)
Box.test(resid.arima01.namer.resid, lag = 10, type = c("Box-Pierce", "Ljung-Box"), fitdf = 3)

# plot the influx data and combined fitted values
influx.fit.namer <- influx.gam01.namer.fit+resid.arima01.namer.fit
plot(influx.ts.namer, type="o")
lines(influx.fit.namer, col="purple", lwd=2)

# poor man's R^2
R <- cor(influx.ts.namer, influx.gam01.namer.fit+resid.arima01.namer.fit)
Rsq <- R^2
print (c(R,Rsq))

# write out some data
out.fit <- cbind(fts.namer$agets,tann.ts.namer,influx.ts.namer,influx.gam01.namer.fit,influx.gam01.namer.resid,
                 resid.arima01.namer.fit,resid.arima01.namer.resid,influx.fit.namer) 
write.table(out.fit, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")

# extract components of GAM plot
fv <- predict(influx.gam01.namer,type="terms", se.fit=TRUE) ## get term estimates       
partial.resid.namer <- residuals(influx.gam01.namer,type="working") + fv$fit #fv[,1]
se.fit.namer <- fv$se.fit
se.upper.namer <- fv$fit + 2.0*fv$se.fit
se.lower.namer <- fv$fit - 2.0*fv$se.fit

# check extraction by overplotting
plot(influx.gam01.namer, ylim=c(-0.8,0.4), residuals=TRUE, pch=16 ) # GAM plot
lines(fts.namer$tann, fv$fit, col="red", lwd=0.2)
lines(fts.namer$tann, se.upper.namer, col="red", lwd=0.2)
lines(fts.namer$tann, se.lower.namer, col="red", lwd=0.2)
points(tann.ts.namer, partial.resid.namer, pch=16, col="red", cex=0.4) # partial residuals

# write GAM plot components
gamplt <- cbind(fts.namer$tann,fv$fit,se.upper.namer,se.lower.namer,partial.resid.namer)          
colnames(gamplt) <- c("namer.tann","namer.tann_fit","namer.tann.upper",
                           "namer.tann.lower","namer.tann.partialresid")
write.table(gamplt, file=gamplotfile, row.names=FALSE, col.names=TRUE, sep=",")