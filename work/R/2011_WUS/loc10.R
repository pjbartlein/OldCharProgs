impd <- read.csv("e:\\Projects\\IMPD\\2010\\data\\prop_scar_ann_w104.csv")
impd <- impd[impd$Year >= 1468,]
impd <- impd[impd$Year <= 1990,]

hist(impd$PropScar, breaks=20)
plot(impd$Year, impd$PropScar)
plot(impd$Year, impd$NScar)
plot(impd$Year, impd$NReport)
summary(impd)

newYear <- data.frame(Year=seq(1000, 2000,  1))

# locfit
# library(locfit)

nn_param <- .3

nrep_loc01 <- locfit(NReport ~ lp(Year, deg=1, nn=nn_param),  family="qpoisson", data=impd)
nrep_pred01 <- predict(nrep_loc01, newdata=newYear, se.fit=T)
summary(nrep_loc01)
plot(impd$NReport ~ impd$Year, col="gray")
lines(newYear$Year, nrep_pred01$fit, col="blue", lwd=3)
lines(newYear$Year, nrep_pred01$fit + 2.0*nrep_pred01$se.fit, col="light blue")
lines(newYear$Year, nrep_pred01$fit - 2.0*nrep_pred01$se.fit, col="light blue")

nscar_loc01 <- locfit(NScar ~ lp(Year, deg=1, nn=nn_param),  family="qpoisson", data=impd)
nscar_pred01 <- predict(nscar_loc01, newdata=newYear, se.fit=T)
summary(nscar_loc01)
plot(impd$NScar ~ impd$Year, col="gray")
lines(newYear$Year, nscar_pred01$fit, col="red", lwd=3)
lines(newYear$Year, nscar_pred01$fit + 2.0*nscar_pred01$se.fit, col="pink")
lines(newYear$Year, nscar_pred01$fit - 2.0*nscar_pred01$se.fit, col="pink")

propscar_loc01 <- locfit(PropScar ~ lp(Year, deg=2, nn_param <- .3, acri="GCV"),
    family="binomial", data=impd)
summary(propscar_loc01)
propscar_pred01 <- predict(propscar_loc01, newdata=newYear, se.fit=T)
plot(impd$PropScar ~ impd$Year, col="gray")
lines(newYear$Year, propscar_pred01$fit + 2.0*propscar_pred01$se.fit, col="plum")
lines(newYear$Year, propscar_pred01$fit - 2.0*propscar_pred01$se.fit, col="plum")
lines(newYear$Year, propscar_pred01$fit, col="purple", lwd=3)
