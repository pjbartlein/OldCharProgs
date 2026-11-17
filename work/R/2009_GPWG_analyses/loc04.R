library(pscl)
#attach(impd)

summary(impd)


hurdle(formula, data, subset, na.action, weights, offset,
dist = c("poisson", "negbin", "geometric"),
zero.dist = c("binomial", "poisson", "negbin", "geometric"),
link = c("logit", "probit", "cloglog", "cauchit", "log"),
control = hurdle.control(...),
model = TRUE, y = TRUE, x = FALSE, ...)

h1 <- hurdle(NScar ~ Year, weights=NReport, zero.dist="binomial", link="logit")

h1 <- hurdle(NScar ~ Year,  zero.dist="binomial")

