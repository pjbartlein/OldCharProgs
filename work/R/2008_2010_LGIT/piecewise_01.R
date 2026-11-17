# piecewise regression fits
library(splines)

year <- seq(10000,15000, by=50)

sm_influx <- read.csv("e:/Projects/igbpfire/yd/data/files/yd_influx_zts_100_boot_all_sites.csv")
attach(sm_influx)
plot(yd_influx_zts_100 ~ targ_age, type="l")


sm_influx.bs1 <- lm(yd_influx_zts_100 ~ bs(targ_age, knots=c(14700., 12900., 11700), degree=1))
plot(yd_influx_zts_100 ~ targ_age, type="l")
lines(predict(sm_influx.bs1, data.frame(targ_age=year)) ~ year, col="red")

#plot(influx.bs1$residuals ~ influx.bs1$fitted.values)
#plot(influx.bs1$residuals ~ targ_age)

sm_influx.bs2 <- lm(yd_influx_zts_100 ~ bs(targ_age, knots=c( 12875., 11660), degree=1))
plot(yd_influx_zts_100 ~ targ_age, type="l")
lines(predict(sm_influx.bs2, data.frame(targ_age=year)) ~ year, col="blue")

detach(sm_influx)

influx <- read.csv("e:/Projects/igbpfire/yd/data/files/yd_influx_zts_b.csv")
attach(influx)
plot(ztinflux ~ est_age)

influx.bs1 <- lm(ztinflux ~ bs(est_age, knots=c(14700., 12900., 11700), degree=1))
plot(ztinflux ~ est_age)
lines(predict(influx.bs1, data.frame(est_age=year)) ~ year, col="red")

influx.bs2 <- lm(ztinflux ~ bs(est_age, knots=c(12875., 11660), degree=1))
plot(ztinflux ~ est_age)
lines(predict(influx.bs2, data.frame(est_age=year)) ~ year, col="red")



plot(predict(sm_influx.bs1, data.frame(targ_age=year)), predict(influx.bs1, data.frame(est_age=year))

plot(yd_influx_zts_100 ~ targ_age, type="l")
lines(predict(sm_influx.bs1, data.frame(targ_age=year)) ~ year, col="red")
lines(predict(influx.bs1, data.frame(est_age=year)) ~ year, col="purple")

cbind(year,predict(sm_influx.bs1, data.frame(targ_age=year)))

sm_influx_fit <- data.frame(cbind(year,predict(sm_influx.bs1, data.frame(targ_age=year))))
names(sm_influx_fit) <- c("year", "influx_fit")
write.table(sm_influx_fit, file="e:/Projects/igbpfire/yd/data/files/yd_influx_ztz_100_boot_pfit.csv", 
	row.names=FALSE, sep=",")
	
sm_influx_fit2 <- data.frame(cbind(year,predict(sm_influx.bs2, data.frame(targ_age=year)), 
	predict(influx.bs2, data.frame(est_age=year)) ))
names(sm_influx_fit2) <- c("year", "sm_influx_fit", "influx_fit")
write.table(sm_influx_fit2, file="e:/Projects/igbpfire/yd/data/files/yd_influx_ztz_100_boot_pfit2.csv", 
row.names=FALSE, sep=",")
	
	


sm_peaks <- read.csv("e:/Projects/igbpfire/yd/data/files/yd_peaks_density_100_512_boot_all_peaks_b.csv")
attach(sm_peaks)
plot(yd_peaks_density_100 ~ targ_age, type="l")

sm_peaks.bs1 <- lm(yd_peaks_density_100 ~ bs(targ_age, knots=c(14700., 12900., 11700), degree=1))
plot(yd_peaks_density_100 ~ targ_age, type="l")
lines(predict(sm_peaks.bs1, data.frame(targ_age=year)) ~ year, col="blue")

sm_peaks.bs2 <- lm(yd_peaks_density_100 ~ bs(targ_age, knots=c(12900., 11700), degree=1))
plot(yd_peaks_density_100 ~ targ_age, type="l")
lines(predict(sm_peaks.bs2, data.frame(targ_age=year)) ~ year, col="blue")

sm_peaks.bs3 <- lm(yd_peaks_density_100 ~ bs(targ_age, knots=c(12900.), degree=1))
plot(yd_peaks_density_100 ~ targ_age, type="l")
lines(predict(sm_peaks.bs3, data.frame(targ_age=year)) ~ year, col="blue")

plot(yd_peaks_density_100 ~ targ_age, type="l")
lines(predict(sm_peaks.bs2, data.frame(targ_age=year)) ~ year, col="red")
lines(predict(sm_peaks.bs3, data.frame(targ_age=year)) ~ year, col="purple")


sm_peaks_fit2 <- data.frame(cbind(year,predict(sm_peaks.bs2, data.frame(targ_age=year))))
names(sm_peaks_fit2) <- c("year", "density_fit")
write.table(sm_peaks_fit2, file="e:/Projects/igbpfire/yd/data/files/yd_density_100_512_boot_pfit2.csv", 
	row.names=FALSE, sep=",")
	
sm_peaks_fit3 <- data.frame(cbind(year,predict(sm_peaks.bs3, data.frame(targ_age=year))))
names(sm_peaks_fit3) <- c("year", "density_fit")
write.table(sm_peaks_fit3, file="e:/Projects/igbpfire/yd/data/files/yd_density_100_512_boot_pfit3.csv", 
row.names=FALSE, sep=",")



sm_poll <- read.csv("e:/Projects/igbpfire/yd/data/files/yd_poll_zts_200_boot_all_sites.csv")
attach(sm_poll)
plot(yd_pol_zts_200 ~ targ_age, type="l")

sm_poll.bs1 <- lm(yd_pol_zts_200 ~ bs(targ_age, knots=c(14700., 12900., 11700), degree=1))
plot(yd_pol_zts_200 ~ targ_age, type="l")
lines(predict(sm_poll.bs1, data.frame(targ_age=year)) ~ year, col="green")

sm_poll_fit <- data.frame(cbind(year,predict(sm_poll.bs1, data.frame(targ_age=year))))
names(sm_poll_fit) <- c("year", "pollen_fit")
write.table(sm_poll_fit, file="e:/Projects/igbpfire/yd/data/files/yd_poll_zts_200_boot_pfit1.csv", 
row.names=FALSE, sep=",")



attach(influx)
# create predictor variables for piecewise regression
x1 <- est_age
x21 <- rep(0,length(est_age))
x22 <- rep(0,length(est_age))
x31 <- rep(0,length(est_age))
x32 <- rep(0,length(est_age))
x21[which(x1 >= 12875.0)] <- 1
x22[which(x1 >= 12875.0)] <- est_age[which(x1 >= 12875.0)]
x31[which(x1 >= 11660 & x1 < 12875.0)] <- 1
x32[which(x1 >= 11660 & x1 < 12875.0)] <- est_age[which(x1 >= 11660 & x1 < 12875.0)]
data1 <- data.frame(cbind(ztinflux, x1, x21, x22, x31, x32))

# and for the fitted values line
x1 <- year
x21 <- rep(0,length(year))
x22 <- rep(0,length(year))
x31 <- rep(0,length(year))
x32 <- rep(0,length(year))
x21[which(x1 >= 12875.0)] <- 1
x22[which(x1 >= 12875.0)] <- year[which(x1 >= 12875.0)]
x31[which(x1 >= 11660 & x1 < 12875.0)] <- 1
x32[which(x1 >= 11660 & x1 < 12875.0)] <- year[which(x1 >= 11660 & x1 < 12875.0)]
datap <- data.frame(cbind(x1, x21, x22, x31, x32))

influx.lm3 <- lm(ztinflux ~ x1+x21+x22+x31+x32, data=data1)
summary(influx.lm3)

plot(ztinflux ~ est_age)
lines(predict(influx.lm3, datap) ~ year, col="red")

influx.lm4 <- lm(ztinflux ~ x1+x21+x31+x32, data=data1)
summary(influx.lm4)

#plot(ztinflux ~ est_age)
lines(predict(influx.lm4, datap) ~ year, col="blue")

influx.bs6 <- lm(ztinflux ~ bs(est_age, knots=c(12875., 11660), degree=2))
summary(influx.bs6)
plot(ztinflux ~ est_age)
lines(predict(influx.bs2, data.frame(est_age=year)) ~ year, col="blue")

# redo the dummy variable regressions with each interval in turn as the base regression

# base regression is est_age < 11660
# create predictor variables for piecewise regression
x1 <- est_age
x21 <- rep(0,length(est_age))
x22 <- rep(0,length(est_age))
x31 <- rep(0,length(est_age))
x32 <- rep(0,length(est_age))
x21[which(est_age >= 11660 & est_age < 12875.0)] <- 1
x22[which(est_age >= 11660 & est_age < 12875.0)] <- est_age[which(est_age >= 11660 & est_age < 12875.0)]
x31[which(est_age >= 12875.0)] <- 1
x32[which(est_age >= 12875.0)] <- est_age[which(est_age >= 12875.0)]
data1 <- data.frame(cbind(ztinflux, x1, x21, x22, x31, x32))

# and for the fitted values line
x1 <- year
x21 <- rep(0,length(year))
x22 <- rep(0,length(year))
x31 <- rep(0,length(year))
x32 <- rep(0,length(year))
x21[which(year >= 11660 & year < 12875.0)] <- 1
x22[which(year >= 11660 & year < 12875.0)] <- year[which(year >= 11660 & year < 12875.0)]
x31[which(year >= 12875.0)] <- 1
x32[which(year >= 12875.0)] <- year[which(year >= 12875.0)]
datap <- data.frame(cbind(x1, x21, x22, x31, x32))

influx.lm31 <- lm(ztinflux ~ x1+x21+x22+x31+x32, data=data1)
summary(influx.lm31)
influx.lm31$coeff

plot(ztinflux ~ est_age)
lines(predict(influx.lm31, datap) ~ year, col="red")

# base regression is est_age >= 11660 & est_age < 12875.0
# create predictor variables for piecewise regression
x11 <- rep(0,length(est_age))
x12 <- rep(0,length(est_age)) 
x2 <- est_age
x31 <- rep(0,length(est_age))
x32 <- rep(0,length(est_age))
x11[which(est_age <= 11660.0)] <- 1
x12[which(est_age <= 11660.0)] <- est_age[which(est_age < 11660.0)]
x31[which(est_age > 12875.0)] <- 1
x32[which(est_age > 12875.0)] <- est_age[which(est_age > 12875.0)]
data1 <- data.frame(cbind(ztinflux, x11, x12, x2, x31, x32))

# and for the fitted values line
x11 <- rep(0,length(year))
x12 <- rep(0,length(year))
x2 <- year
x31 <- rep(0,length(year))
x32 <- rep(0,length(year))
x11[which(year <= 11660.0)] <- 1
x12[which(year <= 11660.0)] <- year[which(year <= 11660.0)]
x31[which(year > 12875.0)] <- 1
x32[which(year > 12875.0)] <- year[which(year > 12875.0)]
datap <- data.frame(cbind(x11, x12, x2, x31, x32))

influx.lm32 <- lm(ztinflux ~ x11+x12+x2+x31+x32, data=data1)
summary(influx.lm32)
influx.lm32$coeff

#plot(ztinflux ~ est_age)
lines(predict(influx.lm32, datap) ~ year, col="blue")

# base regression is est_age <= 12875.0
# create predictor variables for piecewise regression
x11 <- rep(0,length(est_age))
x12 <- rep(0,length(est_age)) 
x21 <- rep(0,length(est_age))
x22 <- rep(0,length(est_age))
x3 <- est_age
x11[which(est_age <= 11660.0)] <- 1
x12[which(est_age <= 11660.0)] <- est_age[which(est_age < 11660.0)]
x21[which(est_age >= 11660 & est_age < 12875.0)] <- 1
x22[which(est_age >= 11660 & est_age < 12875.0)] <- est_age[which(est_age >= 11660 & est_age < 12875.0)]
x3 <- est_age
data1 <- data.frame(cbind(ztinflux, x11, x12, x21, x22, x3))

# and for the fitted values line
x11 <- rep(0,length(year))
x12 <- rep(0,length(year))
x21 <- rep(0,length(year))
x22 <- rep(0,length(year))
x32 <- year
x11[which(year <= 11660.0)] <- 1
x12[which(year <= 11660.0)] <- year[which(year <= 11660.0)]
x21[which(year >= 11660 & year < 12875.0)] <- 1
x22[which(year >= 11660 & year < 12875.0)] <- year[which(year >= 11660 & year < 12875.0)]
x3 <- year
datap <- data.frame(cbind(x11, x12, x21, x22, x3))

influx.lm33 <- lm(ztinflux ~ x11+x12+x21+x22+x3, data=data1)
summary(influx.lm33)
influx.lm33$coeff

#plot(ztinflux ~ est_age)
lines(predict(influx.lm33, datap) ~ year, col="purple")
