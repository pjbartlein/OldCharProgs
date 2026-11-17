# wus influx  -- linear (polynomial) models

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

lm.null.wus <- lm(influx ~ 1, data=wus)
summary(lm.null.wus)

lm.wus <- lm(influx ~ tann + dai + I(dai^2) + I(tann*dai), data=wus )

summary(lm.wus)

plot(lm.wus)
r <- cor(wus$influx, lm.wus$fitted.values); r2 <- r^2; r; r2
anova(lm.null.wus, lm.wus,test="F")

new.tann <- seq(from=-0.3, to=0.1, len=40); new.pdsi <- seq(from=-0.3, to=0.0, len=38)
new.dai <- seq(from=0.3, to=0.4, len=40)
##new.wus$tmp <- seq(from=min(wus$tmp), to=max(wus$tmp), len=40); new.wus$pme <- seq(from=min(wus$pme), to=max(wus$pme), len=40)
new.x <- expand.grid(tann=new.tann, dai=new.dai)
lm.sfc <- predict(lm.wus, new.x)
lm.sfc <- matrix(lm.sfc, 40, 40)
par(mfcol=c(1,1))
persp(new.tann, new.dai, lm.sfc, theta=-60, phi=20, d=1.5, col="gray",ticktype="detailed",
      zlim=c(-0.5,0.5),
      xlab="MAT (C)", ylab="DAI (JJA)",
      zlab = "Z-Score of Transformed Charcoal Influx")

wus.fit <- predict(lm.wus)
wusin.fit <- predict(lm.wus, wusin)
plot(wusin$influx ~ wusin$Year, col="red", pch=16, ylim=c(-.5,.5))
points(wus.fit ~ wus$Year, col="purple", pch=16, cex=.7)
points(wusin.fit ~ wusin$Year, col="blue")

out.wus.fit <- cbind(wusin$Year,wusin$influx,wusin.fit)
colnames(out.wus.fit) <- c("Year","Influx","Pred")
#write.table(out.wus.fit, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")