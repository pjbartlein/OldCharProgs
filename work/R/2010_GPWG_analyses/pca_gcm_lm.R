# regressions of PCs on zonal averages
infile <-"e:\\Projects\\GPWG\\analysis\\output\\v2n\\v2n_pcafiles\\scores_gcm.csv"
outfile <-"e:\\Projects\\GPWG\\analysis\\output\\v2n\\v2n_pcafiles\\scores_gcm_fitted_values.csv"

pg <- read.csv(infile)

# regressions
lm.pc1 <- lm(PC1 ~ tmp_globe + pme_globe, data=pg)
summary(lm.pc1)
plot(pg$age, pg$PC1, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="PC1")
lines(pg$age, lm.pc1$fitted.values, lw=3, col="purple")

lm.pc2 <- lm(PC2 ~ tmp_globe + pme_globe, data=pg)
summary(lm.pc2)
plot(pg$age, pg$PC2, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="PC2")
lines(pg$age, lm.pc2$fitted.values, lw=3, col="purple")



out <- cbind(pg$age,lm.pc1$fitted.values,lm.pc2$fitted.values)

write.table(out, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")