# local regressions of ice-core data

datapath <- "e:/Projects/GPWG/work/feedback/data/ice_core/"

datafile <- "sapart_dc13_krcorr.csv"
nice <- read.csv(paste(datapath,datafile,sep=""))

summary(nice)

plot (nice$d13C ~ nice$Age_yrBP1950, xlim=c(2000,0))

loess.dC13 <- loess(nice$d13C ~ nice$Age_yrBP1950, span=0.35, order=2)
loess.dC13
hat <- predict(loess.dC13)
lines(nice$Age_yrBP1950[order(nice$Age_yrBP1950)], hat[order(nice$Age_yrBP1950)], col="red")

datafile <- "southern_dc13.csv"
sice <- read.csv(paste(datapath,datafile,sep=""))

summary(sice)

plot (sice$d13C ~ sice$Age_yrBP1950, xlim=c(2000,0))
loess.dC13 <- loess(sice$d13C ~ sice$Age_yrBP1950, span=0.33, order=2)
loess.dC13
hat <- predict(loess.dC13)
lines(sice$Age_yrBP1950[order(sice$Age_yrBP1950)], hat[order(sice$Age_yrBP1950)], col="red")


plot (sice$d13C ~ sice$Age_yrBP1950, xlim=c(2000,0))
ss.d13C.S <- smooth.spline(sice$Age_yrBP1950, sice$d13C, spar=0.995)
ss.d13C.S
hat <- predict(ss.d13C.S)
lines(sice$Age_yrBP1950[order(sice$Age_yrBP1950)], hat$y[order(sice$Age_yrBP1950)], col="blue")