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

plot(impd$Year, NScar.glm$residual, pch=16, cex=.5)

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

plot(impd$Year, impd$PropScar, xlab="Year", ylab="PropScars")
lines(newYear$Year,(NScar.pred$fit+1.96*NScar.pred$se.fit)/(NReport.pred$fit+1.96*NReport.pred$se.fit), lwd=2, col="plum")
lines(newYear$Year,(NScar.pred$fit-1.96*NScar.pred$se.fit)/(NReport.pred$fit-1.96*NReport.pred$se.fit), lwd=2, col="plum")
lines(newYear$Year,NScar.pred$fit/NReport.pred$fit, lwd=2, col="purple")

