library(locfit)
#attach(impd)

impd2 <- impd[impd$Year>=1550,]

wgt <- impd2$NReport
#wgt[impd2$NScar == 0] <- 0.0
#wgt[which(impd2$NScar == 0)] <- PropScar[which(impd2$NScar == 0)]*max(NReport)
wgt[impd2$NScar < 1] <- 0.0
#plot(wgt ~ impd2$Year)

# logit model
lfit6 <- locfit(NScar ~ lp(Year, deg=1, h=100.0, acri="GCV"), data=impd2, weights=wgt, family="qbinomial", link="logit")
plot(impd$PropScar ~ impd$Year, pch=16, cex=0.5)
plot(lfit6, add=T, band="global")
lines(lfit6, col="purple", lwd=3)
pred_lfit6 <- predict.locfit(lfit6, where="data", set.fit=T, band="global")
#preplot_lfit6 <- preplot.locfit(lfit6, where="data", band="global")
#lines(impd2$Year, wgt/max(wgt)/2)


y <- pred_lfit6$fit
tr <- pred_lfit6$trans
xx <- impd2$Year
yy <- tr(y)
crit <- pred_lfit6$critval$crit.val
cup <- tr((y + crit * pred_lfit6$se.fit))
clo <- tr((y - crit * pred_lfit6$se.fit))
plot(impd$PropScar ~ impd$Year, pch=16, cex=0.5)
lines(xx, yy, col="green", lwd=3)
lines(xx, clo, col="green", lwd=1, lty=2)
lines(xx, cup, col="green", lwd=1, lty=2)


lfit <- cbind(impd2$Year, yy, clo, cup)
colnames(lfit) <- c("Year","prop_fit02","prop_lower","prop_upper")
write.table(lfit, file=file.choose(), row.names=F, col.names=T, sep=",")
