# piecewise regression fits
library(splines)

workpath <-"e:\\Projects\\GPWG\\analysis\\output\\v2n\\v2n_curves2\\"
curvesname <- "v2n_lgit_0100_b.csv"

infile <- paste(workpath,curvesname,sep="")
outfile <- paste(workpath,"v2n_lgit_0100_b_pwfits.csv",sep="")

lgit <- read.csv(infile)
lgit <- lgit[lgit$targ_age >= 9000.0 & lgit$targ_age <= 16000.0,]

plot(lgit$na ~ lgit$targ_age, xlim=c(15000,10000), type="p", pch=16)

lm_na <- lm(lgit$na ~ bs(lgit$targ_age, knots=c(14700, 12900, 11700), degree=1))
summary(lm_na)
na_fit <- predict(lm_na, data.frame(lgit$targ_age))
plot(lgit$na ~ lgit$targ_age, xlim=c(15000,10000), type="p", pch=16)
lines(na_fit ~ lgit$targ_age, col="red", lwd=3)

lm_eur <- lm(lgit$eur ~ bs(lgit$targ_age, knots=c(14700, 12900, 11700), degree=1))
summary(lm_na)
eur_fit <- predict(lm_eur, data.frame(lgit$targ_age))
plot(lgit$eur ~ lgit$targ_age, xlim=c(15000,10000), type="p", pch=16)
lines(eur_fit ~ lgit$targ_age, col="red", lwd=3)

lm_ssa <- lm(lgit$ssa ~ bs(lgit$targ_age, knots=c(14200, 11800), degree=1))
summary(lm_na)
ssa_fit <- predict(lm_ssa, data.frame(lgit$targ_age))
plot(lgit$ssa ~ lgit$targ_age, xlim=c(15000,10000), type="p", pch=16)
lines(ssa_fit ~ lgit$targ_age, col="red", lwd=3)

lm_smainaustl <- lm(lgit$smainaustl ~ bs(lgit$targ_age, knots=c(14200, 12800, 11800), degree=1))
summary(lm_na)
smainaustl_fit <- predict(lm_smainaustl, data.frame(lgit$targ_age))
plot(lgit$smainaustl ~ lgit$targ_age, xlim=c(15000,10000), type="p", pch=16)
lines(smainaustl_fit ~ lgit$targ_age, col="red", lwd=3)

lm_na_1 <- lm(lgit$na ~ bs(lgit$targ_age, knots=c(14700, 12900, 11700), degree=1))
summary(lm_na_1)
na_fit_1 <- predict(lm_na_1, data.frame(lgit$targ_age))
lm_na_2 <- lm(lgit$na ~ bs(lgit$targ_age, knots=c(14700,  11700), degree=1))
summary(lm_na_2)
na_fit_2 <- predict(lm_na_2, data.frame(lgit$targ_age))
AIC(lm_na_1,lm_na_2)
plot(lgit$na ~ lgit$targ_age, xlim=c(15000,10000), type="p", pch=16)
lines(na_fit_1 ~ lgit$targ_age, col="red", lwd=3)
lines(na_fit_2 ~ lgit$targ_age, col="blue", lwd=3)

lm_ssa_1 <- lm(lgit$ssa ~ bs(lgit$targ_age, knots=c(14700, 12900, 11700), degree=1))
summary(lm_ssa_1)
ssa_fit_1 <- predict(lm_ssa_1, data.frame(lgit$targ_age))
lm_ssa_2 <- lm(lgit$ssa ~ bs(lgit$targ_age, knots=c(14700,  11700), degree=1))
summary(lm_ssa_2)
ssa_fit_2 <- predict(lm_ssa_2, data.frame(lgit$targ_age))
AIC(lm_ssa_1,lm_ssa_2)
plot(lgit$ssa ~ lgit$targ_age, xlim=c(15000,10000), type="p", pch=16)
lines(ssa_fit_1 ~ lgit$targ_age, col="red", lwd=3)
lines(ssa_fit_2 ~ lgit$targ_age, col="blue", lwd=3)

lm_smainaustl_1 <- lm(lgit$smainaustl ~ bs(lgit$targ_age, knots=c(14700, 12900, 11700), degree=1))
summary(lm_smainaustl_1)
smainaustl_fit_1 <- predict(lm_smainaustl_1, data.frame(lgit$targ_age))
lm_smainaustl_2 <- lm(lgit$smainaustl ~ bs(lgit$targ_age, knots=c(14700,  11700), degree=1))
summary(lm_smainaustl_2)
smainaustl_fit_2 <- predict(lm_smainaustl_2, data.frame(lgit$targ_age))
AIC(lm_smainaustl_1,lm_smainaustl_2)
plot(lgit$smainaustl ~ lgit$targ_age, xlim=c(15000,10000), type="p", pch=16)
lines(smainaustl_fit_1 ~ lgit$targ_age, col="red", lwd=3)
lines(smainaustl_fit_2 ~ lgit$targ_age, col="blue", lwd=3)


