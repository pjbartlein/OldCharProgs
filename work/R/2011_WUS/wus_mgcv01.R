# wus influx
# library(mgcv)

workpath <-"e:\\Projects\\GPWG\\work\\wus\\R\\"
curvesname <- "boot100.csv"
outname <- "boot100_fitted.csv"
outfile <- paste(workpath, outname, sep="")
                 
wusin <- read.csv(paste(workpath, curvesname, sep=""))
summary(wusin)
wus <- wusin[wusin$Year <= 1850,]
#wus <- wusin
summary(wus)

par(mfcol=c(1,1))

# regressions

mgcv.null.wus <- gam(influx ~ 1, data=wus)
summary(mgcv.null.wus)

#mgcv.wus <- gam(influx ~ te(tann, bs="tp", k=2) + te(dai, bs="tp", k=2) + te(tann, dai, bs="tp", k=4) , data=wus )
mgcv.wus <- gam(influx ~ te(tann, bs="tp", k=2) + te(dai, bs="tp", k=4) + te(tann, dai, bs="tp", k=2) , data=wus )

#mgcv.wus <- gam(influx ~ te(tann, bs="tp", k=3) + te(dai, bs="tp", k=3) + te(tann, dai, bs="tp", k=6), data=wus )


summary(mgcv.wus)

plot.gam(mgcv.wus, residuals=T, se=T, shade=T, pch=16, cex=0.5)
vis.gam(mgcv.wus, plot.type="persp", type="response", color="gray", theta=-60)
gam.check(mgcv.wus)
r <- cor(wus$influx, mgcv.wus$fitted.values); r2 <- r^2; r; r2
anova(mgcv.null.wus, mgcv.wus,test="F")

new.tann <- seq(from=-0.3, to=0.1, len=40); new.pdsi <- seq(from=-0.3, to=0.0, len=38)
new.dai <- seq(from=0.3, to=0.40, len=40)
##new.wus$tmp <- seq(from=min(wus$tmp), to=max(wus$tmp), len=40); new.wus$pme <- seq(from=min(wus$pme), to=max(wus$pme), len=40)
new.x <- expand.grid(tann=new.tann, dai=new.dai)
mgcv.sfc <- predict(mgcv.wus, new.x)
#mgcv.sfc[new.x$dai >= 0.39] <- NA
#mgcv.sfc[new.x$dai >= 0.375 & new.x$tann <= -0.05] <- NA
#mgcv.sfc[new.x$dai >= 0.365 & new.x$tann <= -0.1] <- NA
#mgcv.sfc[new.x$dai <= 0.31 & new.x$tann >= -0.075] <- NA
mgcv.sfc <- matrix(mgcv.sfc, 40, 40)

par(mfcol=c(1,1))
persp(new.tann, new.dai, mgcv.sfc, theta=-50, phi=30, d=1.5, col="gray",ticktype="detailed",
      zlim=c(-0.5,0.5),
      xlab="MAT (C)", ylab="DAI (JJA)",
      zlab = "Z-Score of Transformed Charcoal Influx")

image(new.tann, new.dai, mgcv.sfc)
points(wus$tann, wus$dai)

wus.fit <- predict(mgcv.wus, wus)
wusin.fit <- predict(mgcv.wus, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")

out.wus.fit <- cbind(wusin$Year,wusin$influx,wusin.fit)
colnames(out.wus.fit) <- c("Year","Influx","Pred")
write.table(out.wus.fit, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")