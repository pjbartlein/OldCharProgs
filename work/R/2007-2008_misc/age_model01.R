#PREP - READ IN THE FILES
siteage <- read.csv(file.choose()) # age,depth
sitedepth <- read.csv(file.choose()) # depth
sitename <- "Cygnet L."  # site name
coredate <- -50
attach(siteage)

plot(age ~ depth, main=sitenamesite)

#SMOOTHING SPLINE
plot(age ~ depth, main=sitename)
age.spline1 <- smooth.spline((age - coredate) ~ depth, spar=.30)
plot(site(age - coredate) ~ depth)
hat <- predict(age.spline1)
lines(spline(hat$x, hat$y))
targetage <- predict(age.spline1,depth$depth)
write.table (targetage, file.choose(), row.names=FALSE, sep=",", col.names=c("depth", "intage"))

#LINEAR INTERPOLATION
plot(age ~ depth)
hat <- approx(depth, age, xout=dep$depth)
lines(depth$depth, hat$y)
write.table (hat, file.choose(), row.names=FALSE, sep=",", col.names=c("depth", "intage"))

#POLYNOMIAL WITHOUT INTERCEPT
plot(age ~ depth)
age.lm1 <- lm(age ~ -1 + depth + I(depth^2))
summary(age.lm1)
lines(spline(fitted(age.lm1) ~ depth), col="blue")
targetages.lm1 <- predict(age.lm1, newdata=depth)
write.table (cbind(depth$depth,targetages.lm1), file.choose(), row.names=FALSE, sep=",", col.names=c("depth", "intage"))

#POLYNOMIAL WITH INTERCEPT
plot(age ~ depth)
age.lm1 <- lm((age - coredate) ~ depth + I(depth^2) + I(depth^3))
summary(age.lm1)
lines(spline(fitted(age.lm1) ~ depth), col="red")
targetages.lm1 <- predict(age.lm1, newdata=depth)
write.table (cbind(depth$depth,targetages.lm1), file.choose(), row.names=FALSE, sep=",", col.names=c("depth", "intage"))

#CLEANUP
detach(age)

