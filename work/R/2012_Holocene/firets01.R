dataname <- "e:\\Projects\\firets\\data\\euro_influx_tann.csv"
fts <- read.csv(dataname)
attach(fts)

start <- agets[length(agets)]
end <- end=agets[1]
freq <- 1/20
deltat <- 20

#influx.ts <- as.ts(influx, start=start, end=end, deltat=deltat)
#tann.ts <- as.ts(tann, start=start, end=end, deltat=deltat)

influx.ts <- as.ts(influx, start=start, end=end, deltat=deltat)
tann.ts <- as.ts(tann, start=start, end=end, deltat=deltat)

fts.spec <- spectrum(ts.union(influx.ts,tann.ts), kernel("modified.daniell", c(5,7)))
plot(fts.spec)
plot(fts.spec, plot.type = "coherency")
plot(fts.spec, plot.type = "phase")

ar.mod01 <- ar(ts.union(influx.ts,tann.ts), order.max=2)
acf(ar.mod01$resid[,1], na.action=na.pass)
acf(ar.mod01$resid[,2], na.action=na.pass)

arima.mod01 <- arima(influx.ts, order=c(1,1,2), xreg=tann.ts)
arima.mod01
acf(arima.mod01$resid, na.action=na.pass)