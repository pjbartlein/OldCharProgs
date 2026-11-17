# GAM plus ARIMA models of influx, North America and Europe
library(gam)

datafile <- "e:\\Projects\\GPWG\\work\\QSR\\firets\\data\\asia_influx_d18O_200.csv"
outfile <- "e:\\Projects\\GPWG\\work\\QSR\\firets\\data\\asia_200_out.csv"
gamplotfile <- "e:\\Projects\\GPWG\\work\\QSR\\firets\\data\\asia_200_gamplot.csv"

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
Box.test(influx.ts.asia, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 3)

# fit a GAM to influx data using d18O as a predictor 
# using gam library
influx.gam01.asia <- gam(influx.ts.asia ~ lo(d18O.ts.asia, span=0.85))
summary(influx.gam01.asia)
AIC(influx.gam01.asia)
plot(influx.gam01.asia, residuals=TRUE, se=TRUE, pch=16)

# plot the observed and fitted influx values
influx.gam01.asia.fit <- ts(influx.gam01.asia$fitted.values, start=begage,  end=endage, deltat=delta)
influx.gam01.asia.resid <- ts(influx.gam01.asia$residuals, start=begage,  end=endage, deltat=delta)
plot(influx.ts.asia, type="o")
lines(influx.gam01.asia.fit, col="orange", lwd=2)

# poor man's R^2
R <- cor(influx.ts.asia, influx.gam01.asia.fit)
Rsq <- R^2
print (c(R,Rsq))

# plot the residuals, and get ACF, PACF and Ljung-Box statistic
ts.plot(influx.gam01.asia.resid)
acf(influx.gam01.asia.resid)
pacf(influx.gam01.asia.resid)
Box.test(influx.gam01.asia.resid, lag = 20, type = c("Box-Pierce", "Ljung-Box"), fitdf = 2)


# ARIMA of GAM residuals (influx.gam01.resid)

# ARMA(1,1) first (final for Europe, looks ok here)
resid.arima01.asia <- arima(influx.gam01.asia.resid,  order=c(1,0,1)) 
resid.arima01.asia
AIC(resid.arima01.asia)

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
out.fit <- cbind(fts.asia$agets,d18O.ts.asia,influx.ts.asia,influx.gam01.asia.fit,influx.gam01.asia.resid,
                 resid.arima01.asia.fit,resid.arima01.asia.resid,influx.fit.asia) 
write.table(out.fit, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")

# extract components of GAM plot
fv <- predict(influx.gam01.asia,type="terms", se.fit=TRUE) ## get term estimates       
partial.resid.asia <- residuals(influx.gam01.asia,type="working") + fv$fit #fv[,1]
se.fit.asia <- fv$se.fit
se.upper.asia <- fv$fit + 2.0*fv$se.fit
se.lower.asia <- fv$fit - 2.0*fv$se.fit

# check extraction by overplotting
plot(influx.gam01.asia, ylim=c(-0.8,0.4), residuals=TRUE, pch=16 ) # GAM plot
lines(fts.asia$d18O, fv$fit, col="red", lwd=0.2)
lines(fts.asia$d18O, se.upper.asia, col="red", lwd=0.2)
lines(fts.asia$d18O, se.lower.asia, col="red", lwd=0.2)
points(d18O.ts.asia, partial.resid.asia, pch=16, col="red", cex=0.4) # partial residuals

# write GAM plot components
gamplt <- cbind(fts.asia$d18O,fv$fit,se.upper.asia,se.lower.asia,partial.resid.asia)          
colnames(gamplt) <- c("asia.d18O","asia.d18O_fit","asia.d18O.upper",
                      "asia.d18O.lower","asia.d18O.partialresid")
write.table(gamplt, file=gamplotfile, row.names=FALSE, col.names=TRUE, sep=",")
