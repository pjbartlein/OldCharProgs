impd <- read.csv("e:\\Projects\\IMPD\\2010\\data\\prop_scar_ann_w104.csv")
impd <- impd[impd$Year >= 1468,]
impd <- impd[impd$Year <= 1990,]

hist(impd$PropScar, breaks=20)
plot(impd$Year, impd$PropScar)
plot(impd$Year, impd$NScar)
plot(impd$Year, impd$NReport)
summary(impd)

newYear <- data.frame(Year=seq(1000, 2000,  1))

# GLMs
f1 <- formula(NScar ~ poly(Year, 8))
NScar.glm <- glm(f1, family=poisson, data=impd)
NScar.pred <- predict(NScar.glm, newdata=newYear, type="response", se=TRUE)
plot(impd$Year, impd$NScar, xlab="Year", ylab="NScars")
lines(newYear$Year, NScar.pred$fit+1.96*NScar.pred$se.fit, lwd=2, col="pink")
lines(newYear$Year, NScar.pred$fit-1.96*NScar.pred$se.fit, lwd=2, col="pink")
lines(newYear$Year, NScar.pred$fit, lwd=2, col="red")

plot(newYear$Year, NScar.pred$fit-1.96*NScar.pred$se.fit, pch=16, cex=.5)

f1 <- formula(NReport ~ poly(Year, 6))
NReport.glm <- glm(f1, family=poisson, data=impd)
NReport.pred <- predict( NReport.glm, newdata=newYear, type="response", se=TRUE)
plot(impd$Year, impd$NReport, xlab="Year", ylab="NReport")
lines(newYear$Year, NReport.pred$fit+1.96*NReport.pred$se.fit, lwd=2, col="light blue")
lines(newYear$Year, NReport.pred$fit-1.96*NReport.pred$se.fit, lwd=2, col="light blue")
lines(newYear$Year, NReport.pred$fit, lwd=2, col="blue")

f1 <- formula(cbind(NScar, NReport-NScar)  ~ poly(Year, 8))
PropScar.glm <- glm(f1, family=binomial, data=impd)
PropScar.pred <- predict(PropScar.glm, newdata=newYear, type="response", se=TRUE)
plot(impd$Year, impd$PropScar, xlab="Year", ylab="PropScar")
lines(newYear$Year,PropScar.pred$fit+1.96*PropScar.pred$se.fit, lwd=2, col="plum")
lines(newYear$Year,PropScar.pred$fit-1.96*PropScar.pred$se.fit, lwd=2, col="plum")
lines(newYear$Year,PropScar.pred$fit, lwd=2, col="purple")

plot(impd$Year, impd$PropScar - PropScar.glm$fitted.values)

plot(impd$Year, impd$PropScar, xlab="Year", ylab="NScars")
lines(newYear$Year,NScar.pred$fit/NReport.pred$fit, lwd=2, col="purple")



# ZIP and ZINB models -- these basically didn't work

# library(pscl)
# library(lmtest)

newdata <- data.frame(Year=seq(1000, 2008,  1), NReport=rep(200, 1009))

f1 <- formula(NScar ~ poly(Year, 6) | NReport )
Zip1 <- zeroinfl(f1, dist="poisson", link = "logit", data=impd)
summary(Zip1)
NScar.pred <- predict(Zip1, newdata=newdata, type="response", se=TRUE)
plot(impd$Year, impd$NScar, xlab="Year", ylab="NScars")
#lines(newYear$Year, NScar.pred$fit+1.96*NScar.pred$se.fit, lwd=2, col="pink")
#lines(newYear$Year, NScar.pred$fit-1.96*NScar.pred$se.fit, lwd=2, col="pink")
lines(newdata$Year, NScar.pred, lwd=2, col="red")


f1 <- (NScar ~ poly(Year, 8)  | NReport )
h1 <- hurdle(f1, dist="negbin", link = "logit", data=impd)
summary(h1)
h1.pred <- predict(h1, newdata=newdata, type="response", se=TRUE)
plot(impd$Year, impd$NScar, xlab="Year", ylab="NScars")
lines(newdata$Year, h1.pred, lwd=2, col="red")

plot(impd$Year, impd$PropScar, xlab="Year", ylab="NScars")
lines(newdata$Year,NScar.pred/NReport.pred$fit, lwd=2, col="purple")

lrtest(Zip1, Nb1)