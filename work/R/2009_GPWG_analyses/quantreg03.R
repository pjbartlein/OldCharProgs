# quantile regression

library(quantreg)
library(splines)

#impd <- read.csv(file.choose())   # read a site_csv file
#attach(impd)


impd2 <- impd[impd$Year>=1300,]

# splines
df <- 6
deg <- 3

#X <- model.matrix(ps2 ~ bs(yr2 , df = df, degree=deg))
#taus <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
taus <- c( 0.5, 0.95, 0.99)
plotclr <- c("red", "purple", "blue")
out.fit <- matrix(nrow=length(impd2$Year), ncol=length(taus))
out.upper <- matrix(nrow=length(impd2$Year), ncol=length(taus))
out.lower <- matrix(nrow=length(impd2$Year), ncol=length(taus))

plot(impd$Year, impd$PropScar, xlab = "Year", ylab = "Proportion", xlim=c(1000,2010), type="n")
points(impd$Year, impd$PropScar, cex = 0.8, col="gray", pch=16)

for (i in 1:length(taus)) {
	print(taus[i])
	fit <- rq(impd2$PropScar ~ bs(impd2$Year , df = df, degree=deg), tau = taus[i]) 
	#out.fit[,i] <- X %*% fit$coef
	out.fit[,i] <- fit$fitted.values
	pred <- predict(fit, newdata=impd2, interval="confidence", level=0.95)
	out.lower[,i] <- pred[,2]
	out.upper[,i] <- pred[,3]
	lines(impd2$Year , out.fit[,i], col=plotclr[i])
	lines(impd2$Year , out.upper[,i], col=plotclr[i], lty=2)
	lines(impd2$Year , out.lower[,i], col=plotclr[i], lty=2)
	}
	
scar.quant <- cbind(impd2$Year, out.fit, out.lower, out.upper)
colnames(scar.quant) <- c("Year","p50","p95","p99","p50_lower","p95_lower","p99_lower",
	"p50_upper","p95_upper","p99_upper")
write.table(scar.quant, file=file.choose(), row.names=F, col.names=T, sep=",")


ps2 <- PropScar[Year >= 1200]
yr2 <- Year[Year >= 1200]
	
# local linear
hs <- c(32, 64)

#X <- model.matrix(ps2 ~ bs(yr2 , df = df, degree=deg))
#taus <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
taus <- c( 0.5, 0.90, 0.95)
#taus <- 0.95
plotclr <- c("red", "purple", "blue")

out.fit <- matrix(nrow=length(yr2), ncol=length(taus))

plot(Year, PropScar, xlab = "Year", ylab = "Proportion", xlim=c(1000,2010), type="n")
points(Year, PropScar, cex = 0.8, col="gray", pch=16)

for (j in 1:length(hs)) {
	for (i in 1:length(taus)) {
		print(c(hs[j], taus[i]))
		fit <- lprq(yr2, ps2, h = hs[j], tau = taus[i]) 
		lines(fit$xx, fit$fv, lty=j, col=plotclr[i])
		}
}


# lowess
hlep
span=0.2
deg <- 10

#X <- model.matrix(ps2 ~ bs(yr2 , df = df, degree=deg))
#taus <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
taus <- c( 0.5, 0.90, 0.95)
plotclr <- c("red", "purple", "blue")
out.fit <- matrix(nrow=length(yr2), ncol=length(taus))

plot(Year, PropScar, xlab = "Year", ylab = "Proportion", xlim=c(1000,2010), type="n")
points(Year, PropScar, cex = 0.8, col="gray", pch=16)

for (i in 1:length(taus)) {
	print(taus[i])
	fit <- rq(ps2 ~ lo(yr2, span=span, degree=deg), tau = taus[i]) 
	out.fit[,i] <- fit$fitted.values
	lines(yr2 , out.fit[,i], col=plotclr[i])
	}	
	
scar.quant <- cbind(yr2,out.fit)
colnames(scar.quant) <- c("Year","p50","p95","p99")
write.table(scar.quant, file=file.choose(), row.names=F, col.names=T, sep=",")

detach(impd)
	
	
# additive fit (qss)

constraint="N"
lambda=1

#X <- model.matrix(ps2 ~ bs(yr2 , df = df, degree=deg))
#taus <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
taus <- c( 0.5, 0.95, 0.99)
plotclr <- c("red", "purple", "blue")
out.fit <- matrix(nrow=length(impd2$Year), ncol=length(taus))
out.upper <- matrix(nrow=length(impd2$Year), ncol=length(taus))
out.lower <- matrix(nrow=length(impd2$Year), ncol=length(taus))

plot(impd$Year, impd$PropScar, xlab = "Year", ylab = "Proportion", xlim=c(1000,2010), type="n")
points(impd$Year, impd$PropScar, cex = 0.8, col="gray", pch=16)

for (i in 1:length(taus)) {
	print(taus[i])
	fit <- rqss(impd$PropScar ~ qss(impd$Year , lambda=20.0), tau = taus[i]) 
	plot(fit)
	#out.fit[,i] <- X %*% fit$coef
	#out.fit[,i] <- fit$fitted.values
	#pred <- predict(fit, newdata=impd2, interval="confidence", level=0.95)
	#out.lower[,i] <- pred[,2]
	#out.upper[,i] <- pred[,3]
	#lines(impd2$Year , out.fit[,i], col=plotclr[i])
	#lines(impd2$Year , out.upper[,i], col=plotclr[i], lty=2)
	#lines(impd2$Year , out.lower[,i], col=plotclr[i], lty=2)
	}
	