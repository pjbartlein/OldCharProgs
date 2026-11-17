# Revised age-model fitting 
# plots have control points (filled circles), interpolating curve, and
# interpolated values (open circles)

#PREP - READ IN THE FILES
siteage <- read.csv(file.choose()) # age,depth
sitedepth <- read.csv(file.choose()) # depth
sitename <- "Cygnet Lake"  # site name
attach(siteage)
topage <- age[1]
topdepth <- depth[1]

# initial plot
plot(age ~ depth, main=sitename)

#SMOOTHING SPLINE
age.spline1 <- smooth.spline(age ~ depth, spar=.30)
age.spline1
hat <- predict(age.spline1)
lines(spline(hat$x, hat$y), col="purple")
interpage <- predict(age.spline1,sitedepth$depth)
plot(interpage, col="purple", pch=1, cex=.5, main=sitename, sub="Smoothing Spline", xlab="depth", ylab="age")
points(depth, age, col="purple", pch=19)
lines(spline(hat$x, hat$y), lwd=2, col="purple")
write.table (interpage, file.choose(), row.names=FALSE, sep=",", col.names=c("depth", "intage"))

#LINEAR INTERPOLATION (note: cannot extrapolate, fixups will be needed outside of age control range))
hat <- approx(depth, age, xout=sitedepth$depth)
plot(hat, col="green", pch=1, cex=.5, main=sitename, sub="Linear Interpolation")
points(depth, age, col="green", pch=19)
lines(sitedepth$depth, hat$y, col="green")
write.table (hat, file.choose(), row.names=FALSE, sep=",", col.names=c("depth", "intage"))

#POLYNOMIAL WITHOUT INTERCEPT -- constrained to pass through youngest age-depth pair
adjage <- age - topage
adjdepth <- depth - topdepth
sitedepth$adjdepth <- sitedepth$depth - topdepth
#age.lm1 <- lm(adjage) ~ -1 + adjdepth + I(adjdepth^2) # 2nd-degree
age.lm1 <- lm(adjage ~ -1 + adjdepth + I(adjdepth^2) + I(adjdepth^3)) # 3rd-degree
summary(age.lm1)
interpages.lm1 <- predict(age.lm1, newdata=sitedepth) + topage
plot(cbind(sitedepth$depth,interpages.lm1), col="blue", pch=1, cex=.5, main=sitename, sub="Constrained Polynomial", xlab="depth", ylab="age")
lines(spline((fitted(age.lm1) + topage) ~ depth), col="blue")
points(depth, age, col="blue", pch=19)
write.table (cbind(sitedepth$depth,interpages.lm1), file.choose(), row.names=FALSE, sep=",", col.names=c("depth", "intage"))

#POLYNOMIAL WITH INTERCEPT -- not constrained to pass through youngest age-depth pair
#age.lm1 <- lm(age  ~ depth + I(depth^2) # 2nd-degree
age.lm1 <- lm(age  ~ depth + I(depth^2) + I(depth^3)) # 3rd-degree
summary(age.lm1)
interpages.lm1 <- predict(age.lm1, newdata=sitedepth)
plot(cbind(sitedepth$depth,interpages.lm1), col="red", pch=1, cex=.5, main=sitename, sub="Polynomial", xlab="depth", ylab="age")
lines(spline(fitted(age.lm1) ~ depth), col="red")
points(depth, age, col="red", pch=19)
write.table (cbind(sitedepth$depth,interpages.lm1), file.choose(), row.names=FALSE, sep=",", col.names=c("depth", "intage"))

#CLEANUP
detach(siteage)

