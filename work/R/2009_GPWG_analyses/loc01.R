impd <- read.csv(file.choose())
# e:\Projects\IMPD\2008\data\prop_scar.csv

impd_dec <- read.csv(file.choose())
# e:\Projects\IMPD\2008\data\prop_scar_dec.csv

impd_anndec < read.csv(file.choose())
# e:\Projects\IMPD\2008\data\prop_scar_anndec.csv

library(locfit)
#attach(impd)

lfit1 <- locfit(NScar ~ lp(Year, deg=2), weights=NReport, family="binomial", link="logit")
plot(impd$PropScar ~ impd$Year)
plot(lfit1, add=T)

lfit2 <- smooth.lf(Year, PropScar, weights=NReport, family="binomial")
plot(impd$PropScar ~ impd$Year)
lines(lfit2)

lfit3 <- locfit(NScar ~ lp(Year, deg=2, nn=0.2), weights=NReport, family="poisson")
plot(impd$PropScar ~ impd$Year)
plot(lfit3, add=T, col="red")

lfit4 <- locfit(NScar ~ lp(Year, deg=2, nn=0.2), weights=NReport, family="binomial", link="logit")
plot(impd$PropScar ~ impd$Year)
lines(lfit4, col="blue", lwd=3)

impd2 <- impd[Year>=1400,]

lfit4 <- locfit(NScar ~ lp(Year, deg=1, nn=0.2), data=impd, weights=NReport, family="binomial", link="logit")
plot(impd$PropScar ~ impd$Year, pch=16, cex=0.5)
lines(lfit4, acol="blue", lwd=3)

lfit5 <- locfit(NScar ~ lp(Year, deg=1, nn=0.2, acri="GCV"), data=impd2, weights=NReport, family="binomial", link="logit")
plot(impd$PropScar ~ impd$Year, pch=16, cex=0.5)
lines(lfit5, col="green", lwd=3)

lfit5 <- locfit(NScar ~ lp(Year, deg=1, h=100, acri="GCV"), data=impd2, weights=NReport, family="qbinomial")
plot(impd$PropScar ~ impd$Year, pch=16, cex=0.5)
lines(lfit5, col="blue", lwd=3)

# logit model
lfit6 <- locfit(NScar ~ lp(Year, deg=1, h=100, acri="GCV"), data=impd2, weights=NReport, family="qbinomial", link="logit")
plot(impd$PropScar ~ impd$Year, pch=16, cex=0.5)
plot(lfit6, add=T, band="local")
lines(lfit6, col="purple", lwd=3)


lfit7 <- locfit(NScar ~ lp(Year, deg=1, h=100, acri="GCV"), data=impd_dec, weights=NReport, family="qbinomial", link="logit")
plot(impd_dec$PropScar ~ impd_dec$Year, pch=16, cex=0.5)
lines(lfit7, col="red", lwd=3)

lfit8 <- locfit(NScar ~ lp(Year, deg=1, h=100, acri="GCV"), data=impd_anndec, weights=NReport, family="qbinomial", link="logit")
plot(impd_anndec$PropScar ~ impd_anndec$Year, pch=16, cex=0.5)
lines(lfit8, col="green", lwd=3)

impd2 <- impd[Year>=1200,]

wgt <- impd2$NReport
#wgt[impd2$NScar == 0] <- 0.0
#wgt[which(impd2$NScar == 0)] <- PropScar[which(impd2$NScar == 0)]*max(NReport)
wgt[impd2$NScar < 1] <- 0.0
plot(wgt ~ impd2$Year)


NRep2 <- NReport[Year >= 1250 & Year <= 1600]
Year2 <- Year[Year >= 1250 & Year <= 1600]
plot(NRep2 ~ Year2)
nrep2_lm01 <- lm(log(NRep2) ~ Year2)
nrep2_pred01 <- exp(nrep2_lm01$fitted.values)

wgt <- impd2$NReport
#wgt[impd2$Year >= 1250 & impd2$Year <= 1600] <- nrep2_pred01
wgt[impd2$NScar < 0] <- 0.0
#plot(wgt ~ impd2$Year)

# logit model
lfit6 <- locfit(NScar ~ lp(Year, deg=1, h=100.0, acri="GCV"), data=impd2, weights=wgt, family="qbinomial", link="logit")
plot(impd$PropScar ~ impd$Year, pch=16, cex=0.5)
plot(lfit6, add=T, band="global")
lines(lfit6, col="purple", lwd=3)
lines(impd2$Year, wgt/max(wgt)/2)


