impd <- read.csv(file.choose())
# e:\Projects\IMPD\2008\data\prop_scar.csv

impd_dec <- read.csv(file.choose())
# e:\Projects\IMPD\2008\data\prop_scar_dec.csv

impd_anndec < read.csv(file.choose())
# e:\Projects\IMPD\2008\data\prop_scar_anndec.csv

library(locfit)
#attach(impd)

summary(impd)

lfit4 <- locfit(NScar ~ lp(Year, deg=1, nn=0.2), data=impd, weights=NReport, family="binomial", link="logit")
plot(impd$PropScar ~ impd$Year, pch=16, cex=0.5)
lines(lfit4, col="blue", lwd=3)

impd2 <- impd[impd$Year>=1550,]

lfit5 <- locfit(NScar ~ lp(Year, deg=1, nn=0.2, acri="GCV"), data=impd2, weights=NReport, family="binomial", link="logit")
plot(impd$PropScar ~ impd$Year, pch=16, cex=0.5)
lines(lfit5, col="green", lwd=3)

attach(impd)

nrep_fit01 <- locfit(NReport ~ lp(Year, deg=1, nn=0.10, acri="GCV"),  family="poisson")
plot(NReport ~ Year, pch=16, cex=0.5)
plot(nrep_fit01, add=T, band="local")
nrep_pred01 <- predict(nrep_fit01, newdata=Year)
lines(Year, nrep_pred01, col="red", lwd=3)

nscar_fit01 <- locfit(NScar ~ lp(Year, deg=1, nn=0.10, acri="GCV"),  family="poisson")
plot(NScar ~ Year, pch=16, cex=0.5)
plot(nscar_fit01, add=T, band="local")
nscar_pred01 <- predict(nscar_fit01, newdata=Year)
lines(Year, nscar_pred01, col="blue", lwd=3)

nprop <- nscar_pred01/nrep_pred01
plot(impd$PropScar ~ impd$Year, pch=16, cex=0.5)
lines(Year, nprop, col="purple", lwd=3)


library(gamlss)

mod1 <- gamlss(NScar ~ lo(Year, span=0.1), family=ZIP)# fits a constant for mu and sigma
fitted(mod1)[1]
fitted(mod1,"sigma")[1]
fitted.plot(mod1, x=Year)


library(pscl)


hurdle(formula, data, subset, na.action, weights, offset,
dist = c("poisson", "negbin", "geometric"),
zero.dist = c("binomial", "poisson", "negbin", "geometric"),
link = c("logit", "probit", "cloglog", "cauchit", "log"),
control = hurdle.control(...),
model = TRUE, y = TRUE, x = FALSE, ...)

h1 <- hurdle(NScar ~ Year +  I(Year^2), weights=NReport, zero.dist="poisson", link="logit")

h1 <- hurdle(NScar ~ Year,  zero.dist="binomial")

library(vglm)

vfit = vglm(NScar ~ Year+I(Year2), zanegbinomial, trace=TRUE)

predict(fit)[1:9,]
