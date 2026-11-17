# regressions of charcoal curves on gcm zonal averages
# library(mgcv)
workpath <-"e:\\Projects\\GPWG\\analysis\\output\\v2n\\v2n_glfiles\\"
curvesname <- "hemis_stack.csv"
outfile <- paste(workpath, "gamfits_", curvesname, sep="")
                 
stack <- read.csv(paste(workpath, curvesname, sep=""))
npts <- length(stack$age)/4

par(mfcol=c(1,1))

# regressions

age <- stack$age
char <- stack$char
x1 <- stack$tmp
x2 <- stack$pre
x3 <- stack$pme
w <- 1.0/(stack$sefit^2)

# final model
mgcv.stack <- gam(char ~ te(x1, bs="tp", k=4) + te( x3, bs="tp", k=8), weights=w)
summary(mgcv.stack)

# diagnostic plots
par(mfcol=c(1,1))
plot.gam(mgcv.stack, residuals=T, se=T, shade=T, pch=16, cex=0.5, ylim=c(-1,1))
vis.gam(mgcv.stack, plot.type="persp", type="response", color="gray", theta=-60)
gam.check(mgcv.stack)
r <- cor(char, mgcv.stack$fitted.values); r2 <- r^2; r; r2

opar <- par(mfcol=c(2,2))
age2 <- age[1:npts]
nextr <- char[1:npts]; nextr_fit <- mgcv.stack$fitted.values[1:npts]
plot(age2, nextr, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="nextr")
points(age2, nextr_fit, pch=16, cex=0.5, col="deeppink")

ntrop <- char[(npts+1):(npts*2)]; ntrop_fit <- mgcv.stack$fitted.values[(npts+1):(npts*2)]
plot(age2, ntrop, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="ntrop")
points(age2, ntrop_fit, pch=16, cex=0.5, col="deeppink")

strop <- char[(2*npts+1):(npts*3)]; strop_fit <- mgcv.stack$fitted.values[(2*npts+1):(npts*3)]
plot(age2, strop, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="strop")
points(age2, strop_fit, pch=16, cex=0.5, col="deeppink")

sextr <- char[(3*npts+1):(npts*4)]; sextr_fit <- mgcv.stack$fitted.values[(3*npts+1):(npts*4)]
plot(age2, sextr, pch=16, col="gray", ylim=c(-1.5,1.5), xlab="Age", ylab="sextr")
points(age2, sextr_fit, pch=16, cex=0.5, col="deeppink")
par(opar)

out <- cbind(age2,nextr,ntrop,strop,sextr,nextr_fit,ntrop_fit,strop_fit,sextr_fit)
write.table(out, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")

# perspective plot
new.x1 <- seq(from=-4, to=2, len=40); new.x3 <- seq(from=-1, to=1, len=40)
#new.x1 <- seq(from=min(x1), to=max(x1), len=40); new.x3 <- seq(from=min(x3), to=max(x3), len=40)
new.x <- expand.grid(x1=new.x1, x3=new.x3)
mgcv.sfc <- predict(mgcv.stack, new.x)
mgcv.sfc <- matrix(mgcv.sfc, 40, 40)
persp(new.x1, new.x3, mgcv.sfc, theta=-60, phi=20, d=1.5, col="gray",ticktype="detailed",
      zlim=c(-1.2,1), xlab="Temperature", ylab="P-E",
      zlab = "Z-Score of Tranformed Charcoal Influx")

# anova
mgcv.stack.null <- gam(char ~ 1)
summary(mgcv.stack.null)
anova(mgcv.stack.null, mgcv.stack, test="F")
