impd <- read.csv(file.choose())
# e:\Projects\IMPD\2008\data\prop_scar.csv

impd_dec <- read.csv(file.choose())
# e:\Projects\IMPD\2008\data\prop_scar_dec.csv

impd_anndec < read.csv(file.choose())
# e:\Projects\IMPD\2008\data\prop_scar_anndec.csv

library(locfit)


nrep_fit01 <- locfit(NReport ~ lp(Year, deg=1, nn=0.05),  family="qpoisson")
plot(NReport ~ Year, pch=16, cex=0.5)
plot(nrep_fit01, add=T, band="global")
nrep_pred01 <- predict(nrep_fit01, newdata=Year, se.fit=T)
lines(Year, nrep_pred01$fit, col="red", lwd=3)
lines(Year, nrep_pred01$fit + 2.0*nrep_pred01$se.fit, col="blue", lty=3)


nscar_fit01 <- locfit(NScar ~ lp(Year, deg=1, nn=0.10, acri="GCV"),  family="poisson")
plot(NScar ~ Year, pch=16, cex=0.5)
plot(nscar_fit01, add=T, band="local")
nscar_pred01 <- predict(nscar_fit01, newdata=Year, se.fit=T)
lines(Year, nscar_pred01$fit, col="blue", lwd=3)
lines(Year, nscar_pred01$fit + 2.0*nscar_pred01$se.fit, col="blue", lty=3)


nprop <- nscar_pred01/nrep_pred01
plot(impd$PropScar ~ impd$Year, pch=16, cex=0.5)
lines(Year, nprop, col="purple", lwd=3)

NScar_bin <- rep(0.0,length(NScar))
NScar_bin[NScar > 0] <- 1.0

NReport_bin <- rep(0.0,length(NReport))
NReport_bin[NReport > 0] <- 1.0

#NScar_bin <- rep(1,length(NScar))
#NScar_bin[NScar > 0] <- 0

nscar_fit02 <- locfit(NScar_bin ~ lp(Year, deg=1, nn=0.05, acri="GCV"),  family="poisson")
plot(NScar ~ Year, pch=16, cex=0.5)
plot(nscar_fit02, add=T, band="local")
nscar_pred02 <- predict(nscar_fit02, newdata=Year, se.fit=T)
#lines(Year, nscar_pred02$fit, col="blue", lwd=3)
plot(nscar_pred02$fit ~ Year)

nrep_fit02 <- locfit(NReport_bin ~ lp(Year, deg=1, nn=0.05, acri="GCV"),  family="poisson")
plot(NReport ~ Year, pch=16, cex=0.5)
plot(nrep_fit02, add=T, band="local")
nrep_pred02 <- predict(nrep_fit02, newdata=Year, se.fit=T)
lines(Year, nrep_pred02$fit, col="red", lwd=3)

nscar_fit03 <- locfit(NScar ~ lp(Year + NScar_bin, deg=1, nn=0.10, acri="GCV"),  family="poisson")
plot(NScar ~ Year, pch=16, cex=0.5)
plot(nscar_fit03, add=T, band="local")
nscar_pred03 <- predict(nscar_fit03, newdata=Year)
lines(Year, nscar_pred03, col="blue", lwd=3)

YearC <- Year[NScar > 0]
NScarC <- NScar[NScar > 0]
plot(NScarC ~ YearC)

nscar_fit04 <- locfit(NScarC ~ lp(YearC, deg=1, nn=0.125, acri="GCV"),  family="poisson")
#plot(NScar ~ Year, pch=16, cex=0.5)
#plot(nscar_fit04, add=T, band="local")
nscar_pred04 <- predict(nscar_fit04, newdata=Year, se.fit=T)
lines(Year, nscar_pred04$fit, col="blue", lwd=3)

nprop <- nscar_pred04$fit/nrep_pred01$fit
plot(PropScar ~ Year, pch=16, cex=0.5)
lines(Year, nprop, col="purple", lwd=3)

nscar_fit05 <- locfit(NScar ~ lp(Year, deg=1, nn=0.125, acri="GCV"),  family="poisson")
plot(NScar ~ Year, pch=16, cex=0.5)
plot(nscar_fit05, add=T, band="local")
nscar_pred05 <- predict(nscar_fit05, newdata=Year, se.fit=T)
lines(Year, nscar_pred05$fit, col="red", lwd=3)

nprop_05 <- nscar_pred05$fit/nrep_pred01$fit
plot(PropScar ~ Year, pch=16, cex=0.5)
lines(Year, nprop_05, col="green", lwd=3)

NScar_censor <- rep(0,length(NScar))
#NScar_censor[nscar_pred02$fit > 0.5] <- 1
NScar_censor[NScar > 1] <- 1

nscar_fit06 <- locfit(NScar ~ lp(Year, deg=1, nn=0.175, acri="GCV"), cen=NScar_censor,  family="poisson")
plot(NScar ~ Year, pch=16, cex=0.5)
plot(nscar_fit06, add=T, band="global")
nscar_pred06 <- predict(nscar_fit06, newdata=Year, se.fit=T)
lines(Year, nscar_pred06$fit, col="blue", lwd=3)

nprop_06 <- nscar_pred06$fit/nrep_pred01$fit
plot(PropScar ~ Year, pch=16, cex=0.5, ylim=c(0,0.3))
lines(Year, nprop_06, col="blue", lwd=3)

nprop_ave <- ((1.0-nscar_pred02$fit)*nprop_05 + nscar_pred02$fit*nprop_06)
plot(PropScar ~ Year, pch=16, cex=0.5)
lines(Year, nprop_ave, col="purple", lwd=3)


nscar_fit07 <- locfit(NScar ~ lp(Year, deg=1, h=100.0, acri="GCV"),  family="poisson")
plot(NScar ~ Year, pch=16, cex=0.5)
plot(nscar_fit07, add=T, band="global")
nscar_pred07 <- predict(nscar_fit07, newdata=Year, se.fit=T)
lines(Year, nscar_pred07$fit, col="blue", lwd=3)

nprop_07 <- nscar_pred07$fit/nrep_pred01$fit
plot(PropScar ~ Year, pch=16, cex=0.5, ylim=c(0,0.3))
lines(Year, nprop_07, col="blue", lwd=3)
