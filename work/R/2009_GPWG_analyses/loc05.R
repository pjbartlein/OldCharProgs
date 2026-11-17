library(mgcv)

attach(impd)

nrep_fit01 <- gam(NReport ~ s(Year), family=quasipoisson()) 
nrep_pred01 <- predict(nrep_fit01, type="response", se.fit=T)
nrep_pred01_upper <- nrep_fit01$fitted.values + 2.0*nrep_pred01$se.fit
nrep_pred01_lower <- nrep_fit01$fitted.values - 2.0*nrep_pred01$se.fit
plot(NReport ~ Year, pch=16, cex=0.5, col="gray")
lines(nrep_fit01$fitted.values ~ Year, col="red", lwd=3)
lines((nrep_pred01_upper) ~ Year, col="red", lwd=1, lty=2)
lines((nrep_pred01_lower) ~ Year, col="red", lwd=1, lty=2)

nscar_fit01 <- gam(NScar ~ s(Year), family=quasipoisson()) 
nscar_pred01 <- predict(nscar_fit01, type="response", se.fit=T)
nscar_pred01_upper <- nscar_fit01$fitted.values + 2.0*nscar_pred01$se.fit
nscar_pred01_lower <- nscar_fit01$fitted.values - 2.0*nscar_pred01$se.fit
plot(NScar ~ Year, pch=16, cex=0.5, col="gray")
lines(nscar_fit01$fitted.values ~ Year, col="red", lwd=3)
lines((nscar_pred01_upper) ~ Year, col="red", lwd=1, lty=2)
lines((nscar_pred01_lower) ~ Year, col="red", lwd=1, lty=2)

nprop01 <- nscar_pred01$fit/nrep_pred01$fit
plot(impd$PropScar ~ impd$Year, pch=16, cex=0.5)
lines(Year, nprop01, col="green", lwd=3)

NScar_bin <- rep(0.0,length(NScar))
NScar_bin[NScar > 0] <- 1.0

NReport_bin <- rep(0.0,length(NReport))
NReport_bin[NReport > 0] <- 1.0


wgt <- (NReport/max(NReport))
#wgt <- (NScar/max(NScar)
wgt1 <- rep(1.0, length(NScar))
#wgt <- NScar_bin
plot(wgt ~ Year)

#NRep2 <- NReport[Year >= 1250 & Year <= 1600]
#Year2 <- Year[Year >= 1250 & Year <= 1600]
#plot(NRep2 ~ Year2)
#nrep2_lm01 <- lm(log(NRep2) ~ Year2)


nrep_fit02 <- gam(NReport ~ s(Year, k=10, bs="tp"), weights=wgt1, family=quasipoisson()) 
nrep_pred02 <- predict(nrep_fit02, type="response", se.fit=T)
nrep_pred02_upper <- nrep_fit02$fitted.values + 2.0*nrep_pred02$se.fit
nrep_pred02_lower <- nrep_fit02$fitted.values - 2.0*nrep_pred02$se.fit
plot(NReport ~ Year, pch=16, cex=0.5, col="gray")
lines(nrep_fit02$fitted.values ~ Year, col="red", lwd=3)
lines((nrep_pred02_upper) ~ Year, col="red", lwd=1, lty=2)
lines((nrep_pred02_lower) ~ Year, col="red", lwd=1, lty=2)

nscar_fit02 <- gam(NScar ~ s(Year), weights=wgt, family=quasipoisson()) 
nscar_pred02 <- predict(nscar_fit02, type="response", se.fit=T)
nscar_pred02_upper <- nscar_fit02$fitted.values + 2.0*nscar_pred02$se.fit
nscar_pred02_lower <- nscar_fit02$fitted.values - 2.0*nscar_pred02$se.fit
plot(NScar ~ Year, pch=16, cex=0.5, col="gray")
lines(nscar_fit02$fitted.values ~ Year, col="red", lwd=3)
lines((nscar_pred02_upper) ~ Year, col="red", lwd=1, lty=2)
lines((nscar_pred02_lower) ~ Year, col="red", lwd=1, lty=2)

nprop02 <- nscar_pred02$fit/nrep_pred02$fit
plot(impd$PropScar ~ impd$Year, pch=16, cex=0.5)
lines(Year, nprop02, col="purple", lwd=3)


