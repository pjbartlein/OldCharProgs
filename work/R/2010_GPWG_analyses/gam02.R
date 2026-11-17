# regressions of charcoal curves on gcm zonal averages
# library(mgcv)
workpath <-"e:\\Projects\\cspace\\data\\work04\\"

cspacename <- "gfed3_mat_pme_globe.csv"
outfile <- paste(workpath, "gamfits_1d_", cspacename, sep="")
                 
csp <- read.csv(paste(workpath, cspacename, sep=""))

csp$gfed3 <- sqrt(csp$gfed3)
#csp$gfed3 <- log10(csp$gfed3+0.001)
#csp$pme <- csp$pme^0.3333

plot(csp$mat, csp$gfed3, pch=16, cex=0.5)
plot(csp$pme, csp$gfed3, pch=16, cex=0.5)

# regressions
mgcv.csp <- gam(gfed3 ~ te(mat, bs="tp", k=4) + te(pme, bs="tp", k=8), data=csp)
summary(mgcv.csp)

# diagnostic plots
par(mfcol=c(1,1))
plot.gam(mgcv.csp, residuals=F, se=T, shade=T, pch=16, cex=0.5, ylim=c(0-.1,.2))
#vis.gam(mgcv.stack, plot.type="persp", type="response", color="gray", theta=-60)
gam.check(mgcv.csp)

# anova
mgcv.csp.null <- gam(sqrt(gfed3) ~ 1, data=csp)
summary(mgcv.csp.null)
anova(mgcv.csp.null, mgcv.csp, test="F")
r <- cor(csp$gfed3, mgcv.csp$fitted.values); r2 <- r^2; r; r2

# perspective plot
par(mfcol=c(1,1))
new.mat <- seq(from=-20, to=30, len=40); new.pme <- seq(from=0.0, to=1500, len=40)
#new.mat <- seq(from=min(csp$mat), to=max(csp$mat), len=40); new.pme <- seq(from=min(csp$pme), to=max(csp$pme), len=40)
new.x <- expand.grid(mat=new.mat, pme=new.pme)
mgcv.sfc <- predict(mgcv.csp, new.x)
mgcv.sfc <- matrix(mgcv.sfc, 40, 40)
persp(new.mat, new.pme, mgcv.sfc, theta=-65, phi=30, d=1.5, col="gray", ticktype="detailed",
      zlim=c(-.1,.3), xlab="MAT (C)", ylab="Annual P-E (mm)", xlim=c(-20,30), ylim=c(0,1500),
      zlab = "GFED 3.1 Burned Fraction")

new.pme2 <- 250.0
new.x <- expand.grid(mat=new.mat, pme=new.pme2)
mgcv.fit.mat <- predict(mgcv.csp, new.x, se.fit="TRUE")
new.mat2 <- 20.0
new.x <- expand.grid(mat=new.mat2, pme=new.pme)
mgcv.fit.pme <- predict(mgcv.csp, new.x, se.fit="TRUE")

mgcv.csp.1dse <- cbind(new.tmp, mgcv.fit.mat$fit, mgcv.fit.mat$fit+2*mgcv.fit.mat$se.fit, mgcv.fit.mat$fit-2*mgcv.fit.mat$se.fit,
                       new.pme, mgcv.fit.pme$fit, mgcv.fit.pme$fit+2*mgcv.fit.pme$se.fit, mgcv.fit.pme$fit-2*mgcv.fit.pme$se.fit) 
colnames(mgcv.csp.1dse) <- c("mat","matfit","matfitlcl","matfitucl","pme","pmefit","pmelcl","pmeucl")
write.table(mgcv.csp.1dse, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")

#stack$ztp <- ((stack$zt-min(stack$zt))/(max(stack$zt)-min(stack$zt))+.01)^0.33
