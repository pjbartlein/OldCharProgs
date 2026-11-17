# wus influx
# library(locfit)

workpath <-"e:\\Projects\\GPWG\\work\\wus\\R\\"
curvesname <- "boot100.csv"
outname <- "boot100_fitted_1800.csv"
outfile <- paste(workpath, outname, sep="")
                 
wusin <- read.csv(paste(workpath, curvesname, sep=""))
summary(wusin)
wus <- wusin[wusin$Year <= 1800,]
#wus <- wusin
summary(wus)

par(mfcol=c(1,1))

# regressions -- 

locfit.wus.null <- lm(influx ~ 1, data=wus)
summary(locfit.wus.null)

# locfit
locfit.wus <- locfit(influx ~ tann + dai, alpha=0.95, data=wus )
summary(locfit.wus)
plot(locfit.wus, get.data=T)
plot(locfit.wus, type="persp")

wus.fit <- predict(locfit.wus, wus)
wusin.fit <- predict(locfit.wus, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")
r <- cor(wus$influx, wus.fit); r2 <- r^2; r; r2



# plots
new.tann <- seq(from=-0.3, to=0.1, len=40); new.pdsi <- seq(from=-0.3, to=0.0, len=38)
new.dai <- seq(from=0.3, to=0.40, len=40)
##new.wus$tmp <- seq(from=min(wus$tmp), to=max(wus$tmp), len=40); new.wus$pme <- seq(from=min(wus$pme), to=max(wus$pme), len=40)
new.x <- expand.grid(tann=new.tann, dai=new.dai)
lm.sfc <- predict(locfit.wus, new.x)
#lm.sfc[new.x$dai >= 0.39] <- NA
#lm.sfc[new.x$dai >= 0.375 & new.x$tann <= -0.05] <- NA
#lm.sfc[new.x$dai >= 0.365 & new.x$tann <= -0.1] <- NA
#lm.sfc[new.x$dai <= 0.31 & new.x$tann >= -0.075] <- NA
lm.sfc <- matrix(lm.sfc, 40, 40)

par(mfcol=c(1,1))
persp(new.tann, new.dai, lm.sfc, theta=-50, phi=30, d=1.5, col="gray",ticktype="detailed",
      zlim=c(-0.5,0.5),
      xlab="MAT (C)", ylab="DAI (JJA)",
      zlab = "Z-Score of Transformed Charcoal Influx")

image(new.tann, new.dai, lm.sfc)
points(wus$tann, wus$dai)

wus.fit <- predict(locfit.wus, wus)
wusin.fit <- predict(locfit.wus, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue"); lines(wusin.fit ~ wusin$Year, col="blue")

out.wus.fit <- cbind(wusin$Year,wusin$influx,wusin.fit)
colnames(out.wus.fit) <- c("Year","Influx","Pred")
write.table(out.wus.fit, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")