library(gam)

#attach(impd)

# this one
impd2 <- impd[impd$Year>=1250,]

wgt <- impd2$NReport
#wgt[impd2$NScar == 0] <- 1
#wgt[which(impd2$NScar < 2)] <- PropScar[which(impd2$NScar < 2)]*max(NReport)

wgt[impd2$NScar < 3] <- 1

pscar_fit01 <- gam(PropScar ~ lo(Year, span=0.2, degree=2), weights=wgt, family=quasibinomial(), data=impd2) 
#pscar_fit01 <- gam(PropScar ~ s(Year, spar=0.75), weights=wgt, family=quasibinomial(), data=impd2) 
#pscar_fit01 <- gam(PropScar ~ s(Year, df=10), weights=wgt, family=quasibinomial(), data=impd2) 
pscar_pred01 <- predict(pscar_fit01, type="response", se.fit=T)
pscar_pred01_upper <- pscar_fit01$fitted.values + 2.0*pscar_pred01$se.fit
pscar_pred01_lower <- pscar_fit01$fitted.values - 2.0*pscar_pred01$se.fit
plot(PropScar ~ Year, pch=16, cex=0.5, col="gray", data=impd)
lines(pscar_fit01$fitted.values ~ impd2$Year, col="green", lwd=3)
lines((pscar_pred01_upper) ~ impd2$Year, col="green", lwd=1, lty=2)
lines((pscar_pred01_lower) ~ impd2$Year, col="green", lwd=1, lty=2)
#lines(impd2$Year, wgt/max(wgt)/2, col="blue")

#wgt <- rep(1.0, length(impd2$NScar))
#wgt <- (impd2$NReport/max(impd$NReport))
wgt2 <- impd2$NReport
#wgt[impd2$NScar == 0] <- 1
#wgt[which(impd2$NScar < 2)] <- PropScar[which(impd2$NScar < 2)]*max(NReport)
wgt2[impd2$NScar < 3] <- 1

nrep_fit02 <- gam(NReport ~ lo(Year, span=0.2, degree=2), weights=wgt, family=quasipoisson(), data=impd2) 
#nrep_fit02 <- gam(NReport ~ s(Year, df=20), weights=wgt, family=quasipoisson(), data=impd2) 
nrep_pred02 <- predict(nrep_fit02, type="response", se.fit=T)
nrep_pred02_upper <- nrep_fit02$fitted.values + 2.0*nrep_pred02$se.fit
nrep_pred02_lower <- nrep_fit02$fitted.values - 2.0*nrep_pred02$se.fit
plot(NReport ~ Year, pch=16, cex=0.5, col="gray", data=impd)
lines(nrep_fit02$fitted.values ~ impd2$Year, col="red", lwd=3)
lines((nrep_pred02_upper) ~ impd2$Year, col="red", lwd=1, lty=2)
lines((nrep_pred02_lower) ~ impd2$Year, col="red", lwd=1, lty=2)

nscar_fit02 <- gam(NScar ~ lo(Year, span=0.4, degree=2), weights=wgt2, family=quasipoisson(), data=impd2) 
#nscar_fit02 <- gam(NScar ~ s(Year, df=8), weights=wgt2, family=quasipoisson(), data=impd2) 
nscar_pred02 <- predict(nscar_fit02, type="response", se.fit=T)
nscar_pred02_upper <- nscar_fit02$fitted.values + 2.0*nscar_pred02$se.fit
nscar_pred02_lower <- nscar_fit02$fitted.values - 2.0*nscar_pred02$se.fit
plot(NScar ~ Year, pch=16, cex=0.5, col="gray", data=impd)
lines(nscar_fit02$fitted.values ~ impd2$Year, col="red", lwd=3)
lines((nscar_pred02_upper) ~ impd2$Year, col="red", lwd=1, lty=2)
lines((nscar_pred02_lower) ~ impd2$Year, col="red", lwd=1, lty=2)

prop02 <- nscar_pred02$fit/nrep_pred02$fit
prop02_upper <- nscar_pred02_upper/nrep_pred02_upper
prop02_lower <- nscar_pred02_lower/nrep_pred02_lower
plot(impd$PropScar ~ impd$Year, pch=16, cex=0.5)
lines(impd2$Year, prop02, col="purple", lwd=3)
lines((prop02_upper) ~ impd2$Year, col="purple", lwd=1, lty=2)
lines((prop02_lower) ~ impd2$Year, col="purple", lwd=1, lty=2)

nrep <- cbind(impd2$Year, nrep_fit02$fitted.values, nrep_pred02_lower, nrep_pred02_upper)
colnames(nrep) <- c("Year","nrep_fit02","nrep_lower","nrep_upper")
write.table(nrep, file=file.choose(), row.names=F, col.names=T, sep=",")

nscar <- cbind(impd2$Year, nscar_fit02$fitted.values, nscar_pred02_lower, nscar_pred02_upper)
colnames(nscar) <- c("Year","nscar_fit02","nscar_lower","nscar_upper")
write.table(nscar, file=file.choose(), row.names=F, col.names=T, sep=",")

prop <- cbind(impd2$Year, prop02, prop02_lower, prop02_upper)
colnames(prop) <- c("Year","prop_fit02","prop_lower","prop_upper")
write.table(prop, file=file.choose(), row.names=F, col.names=T, sep=",")
