# wus influx
# library(mgcv)

workpath <-"e:\\Projects\\GPWG\\work\\wus\\R\\"
curvesname <- "boot100.csv"
outname <- "boot100_fitted_1800.csv"
outfile <- paste(workpath, outname, sep="")
                 
wusin <- read.csv(paste(workpath, curvesname, sep=""))
summary(wusin)
wus <- wusin[wusin$Year <= 1850,]
#wus <- wusin
summary(wus)

par(mfcol=c(1,1))

# regressions -- 

mgcv.wus.null <- lm(influx ~ 1, data=wus)
summary(mgcv.wus.null)

#mgcv.wus <- gam(influx ~ te(tann, bs="tp", k=2) + te(dai, bs="tp", k=2) + te(tann, dai, bs="tp", k=4) , data=wus )

# tann alone
mgcv.wus.t <- gam(influx ~ s(tann, k=4), data=wus )
summary(mgcv.wus.t)
plot.gam(mgcv.wus.t, residuals=T, se=T, shade=T, pch=16, cex=0.5)
gam.check(mgcv.wus.t)
wus.fit <- predict(mgcv.wus.t, wus)
wusin.fit <- predict(mgcv.wus.t, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, mgcv.wus.t$fitted.values); r2 <- r^2; r; r2
anova(mgcv.wus.null, mgcv.wus.t,test="F")
AIC(mgcv.wus.null, mgcv.wus.t)

# pdsi alone
mgcv.wus.p <- gam(influx ~ s(pdsi, k=4), data=wus )
summary(mgcv.wus.p)
plot.gam(mgcv.wus.p, residuals=T, se=T, shade=T, pch=16, cex=0.5)
gam.check(mgcv.wus.p)
wus.fit <- predict(mgcv.wus.p, wus)
wusin.fit <- predict(mgcv.wus.p, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, mgcv.wus.p$fitted.values); r2 <- r^2; r; r2
anova(mgcv.wus.null, mgcv.wus.p,test="F")
AIC(mgcv.wus.null, mgcv.wus.p)

anova(mgcv.wus.t, mgcv.wus.p,test="F")
AIC(mgcv.wus.t, mgcv.wus.p)

# interaction
mgcv.wus.tp <- gam(influx ~ s(tann, pdsi, k=4), data=wus )
summary(mgcv.wus.tp)
plot.gam(mgcv.wus.tp, residuals=T, se=T, shade=T, pch=16, cex=0.5)
gam.check(mgcv.wus.tp)
wus.fit <- predict(mgcv.wus.tp, wus)
wusin.fit <- predict(mgcv.wus.tp, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, mgcv.wus.tp$fitted.values); r2 <- r^2; r; r2
anova(mgcv.wus.null, mgcv.wus.tp,test="F")
AIC(mgcv.wus.null, mgcv.wus.tp)

anova(mgcv.wus.t, mgcv.wus.tp,test="F")
AIC(mgcv.wus.t, mgcv.wus.tp)

# tann + interaction
mgcv.wus.ttp <- gam(influx ~ s(tann, k=2) + s(tann, pdsi, k=4), data=wus )
summary(mgcv.wus.ttp)
plot.gam(mgcv.wus.ttp, residuals=T, se=T, shade=T, pch=16, cex=0.5)
gam.check(mgcv.wus.ttp)
wus.fit <- predict(mgcv.wus.ttp, wus)
wusin.fit <- predict(mgcv.wus.ttp, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, mgcv.wus.ttp$fitted.values); r2 <- r^2; r; r2
anova(mgcv.wus.null, mgcv.wus.ttp,test="F")
AIC(mgcv.wus.null, mgcv.wus.ttp)

anova(mgcv.wus.tp, mgcv.wus.ttp,test="F")
AIC(mgcv.wus.tp, mgcv.wus.ttp)

# pdsi + interaction
mgcv.wus.ptp<- gam(influx ~ s(pdsi, k=2) + te(tann, pdsi, k=4), data=wus )
summary(mgcv.wus.ptp)
plot.gam(mgcv.wus.ptp, residuals=T, se=T, shade=T, pch=16, cex=0.5)
gam.check(mgcv.wus.ptp)
wus.fit <- predict(mgcv.wus.ptp, wus)
wusin.fit <- predict(mgcv.wus.ptp, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, mgcv.wus.ptp$fitted.values); r2 <- r^2; r; r2
anova(mgcv.wus.null, mgcv.wus.ptp,test="F")
AIC(mgcv.wus.null, mgcv.wus.ptp)

anova(mgcv.wus.tp, mgcv.wus.ptp,test="F")
AIC(mgcv.wus.tp, mgcv.wus.ptp)

anova(mgcv.wus.ttp, mgcv.wus.ptp,test="F")
AIC(mgcv.wus.ttp, mgcv.wus.ptp)

# full model
#mgcv.wus <- gam(influx ~ te(tann, k=2) + te(pdsi, k=2) + te(tann, pdsi, k=4), data=wus ) # 2 2 2
#summary(mgcv.wus)

mgcv.wus <- gam(influx ~  s(tann, k=2) + s(pdsi, k=2)  + te(tann, pdsi, k=4), data=wus )  #
summary(mgcv.wus)

plot.gam(mgcv.wus, residuals=T, se=T, shade=T, pch=16, cex=0.5)
vis.gam(mgcv.wus, plot.type="persp", type="response", color="gray", theta=-60)
gam.check(mgcv.wus)
wus.fit <- predict(mgcv.wus, wus)
wusin.fit <- predict(mgcv.wus, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, mgcv.wus$fitted.values); r2 <- r^2; r; r2
anova(mgcv.wus.null, mgcv.wus,test="F")
AIC(mgcv.wus.null, mgcv.wus)

anova(mgcv.wus.ttp, mgcv.wus,test="F")
AIC(mgcv.wus.ttp, mgcv.wus)

# plots
new.tann <- seq(from=-0.3, to=0.1, len=40); new.pdsi <- seq(from=-0.45, to=-0.05, len=40)
#new.dai <- seq(from=0.3, to=0.40, len=40)
##new.wus$tmp <- seq(from=min(wus$tmp), to=max(wus$tmp), len=40); new.wus$pme <- seq(from=min(wus$pme), to=max(wus$pme), len=40)
new.x <- expand.grid(tann=new.tann, pdsi=new.pdsi) # dai=new.dai)
mgcv.sfc <- predict(mgcv.wus, new.x)
#lm.sfc[new.x$dai >= 0.39] <- NA
mgcv.sfc[new.x$pdsi <= -0.25  & new.x$tann <= -0.25] <- NA
mgcv.sfc[new.x$pdsi >= -0.2   & new.x$tann >= -0.075] <- NA
mgcv.sfc[new.x$pdsi <= -0.38  & new.x$tann <=  0.025] <- NA
mgcv.sfc <- matrix(mgcv.sfc, 40, 40)

par(mfcol=c(1,1))
persp(new.tann, new.pdsi, mgcv.sfc, theta=-60, phi=15, d=1.5, col="gray",ticktype="detailed",
      zlim=c(-0.5,0.5),
      xlab="MAT (C)", ylab="PDSI (JJA)", #"DAI (JJA)",
      zlab = "Z-Score of Transformed Charcoal Influx")

image(new.tann, new.pdsi, mgcv.sfc) #new.dai, mgcv.sfc)
points(wus$tann, wus$pdsi) #wus$dai)

wus.fit <- predict(mgcv.wus, wus)
wusin.fit <- predict(mgcv.wus, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")

out.wus.fit <- cbind(wusin$Year,wusin$influx,wusin.fit)
colnames(out.wus.fit) <- c("Year","Influx","Pred")
write.table(out.wus.fit, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")