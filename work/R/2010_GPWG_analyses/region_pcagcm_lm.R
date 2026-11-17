# regressions of charcoal curves on PCs
workpath <-"e:\\Projects\\GPWG\\analysis\\output\\v2n\\v2n_curves2\\"
curvesname <- "v2n_all_nsa_boot_influx_zt21k_sm10_0500.csv"
gcmscores <- "latzones_prcomp_18_scores_2000_lwr.csv"

scoresfile <- paste(workpath,gcmscores,sep="")
curvesfile <- paste(workpath,curvesname,sep="")
outfile <- paste(workpath,"curvepcafits.csv",sep="")

scoresin <- read.csv(scoresfile)
curvesin <- read.csv(curvesfile)
sc <- cbind(scoresin,curvesin[1:211,])

# regressions
lm.nhemis <- lm(nhemis ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7, data=sc)
summary(lm.nhemis)
plot(sc$age, sc$nhemis, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="NHemis")
lines(sc$age, lm.nhemis$fitted.values, lw=3, col="purple")

lm.shemis <- lm(shemis ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7, data=sc)
summary(lm.shemis)
plot(sc$age, sc$shemis, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="SHemis")
lines(sc$age, lm.shemis$fitted.values, lw=3, col="purple")

lm.globe <- lm(globe ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7, data=sc)
summary(lm.globe)
plot(sc$age, sc$globe, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="Globe")
lines(sc$age, lm.globe$fitted.values, lw=3, col="purple")

lm.ntrop <- lm(ntrop ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7, data=sc)
summary(lm.ntrop)
plot(sc$age, sc$ntrop, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="Northern Tropics")
lines(sc$age, lm.ntrop$fitted.values, lw=3, col="purple")

lm.strop <- lm(strop ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7, data=sc)
summary(lm.strop)
plot(sc$age, sc$strop, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="Southern Tropics")
lines(sc$age, lm.strop$fitted.values, lw=3, col="purple")

lm.nextr <- lm(nextr ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7, data=sc)
summary(lm.nextr)
plot(sc$age, sc$nextr, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="Northern Extratropics")
lines(sc$age, lm.nextr$fitted.values, lw=3, col="purple")

lm.sextr <- lm(sextr ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7, data=sc)
summary(lm.sextr)
plot(sc$age, sc$sextr, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="Southern Extratropics")
lines(sc$age, lm.sextr$fitted.values, lw=3, col="purple")

lm.hemis_diff <- lm(hemis_diff ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7, data=sc)
summary(lm.hemis_diff)
plot(sc$age, sc$hemis_diff, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="Hemispheric Difference")
lines(sc$age, lm.hemis_diff$fitted.values, lw=3, col="purple")

out <- cbind(scoresin,lm.nhemis$fitted.values,lm.shemis$fitted.values,lm.globe$fitted.values,
    lm.ntrop$fitted.values,lm.strop$fitted.values,lm.nextr$fitted.values,lm.sextr$fitted.values,
    lm.hemis_diff$fitted.values)

write.table(out, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")