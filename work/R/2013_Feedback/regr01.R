# regression of influx on dc13 and ch4

datapath <- "e:/Projects/GPWG/work/feedback/data/ice_core/"
datafile <- "v3d_icecore_050.csv"

input <- read.csv(paste(datapath,datafile,sep=""))
charice <- input[input$targ_age >= 250.0 & input$targ_age <= 2000,]
summary(charice)

attach(charice)
plot(dc13,nextr)
plot(ch4,nextr)
plot(globe ~ targ_age, pch=16, col="red")

opar <- par(mar=c(5,4,4,5)+0.1) # space for second axis
plot(nextr ~ targ_age, pch=16, col="red")  # first plot
par(new=T)          # second plot is going to get added to first
plot(dc13 ~ targ_age, pch=16, axes=F, ylab="", col="gray")  # second plot
axis(side=4)   # add axis
mtext(side=4,line=3.8,"dc13")   # add label
par(opar)  # restore plot par

# linear regression 
lm01 <- lm(nextr ~ dc13+ch4, data=charice)
summary(lm01)

# regression with dc13 and ch4 product term
lm02 <- lm(nextr ~ I(dc13*ch4) + ch4, data=charice)
summary(lm02)
plot(lm02, which=c(1,2,4))

# GAMs
library(mgcv)

# null model
mgcv.null.gl <- gam(nextr ~ 1, data=charice)
summary(mgcv.null.gl)

# just dc13
mgcv.gl01 <- gam(nextr ~ te(dc13, bs="tp", k=4) , data=charice )
summary(mgcv.gl01)

plot.gam(mgcv.gl01, residuals=T, se=T, shade=T, pch=16, cex=0.5)
gam.check(mgcv.gl01)
anova(mgcv.null.gl, mgcv.gl01,test="F")


# dc13 and ch4
mgcv.gl02 <- gam(globe ~ te(dc13, bs="tp", k=4) + te(ch4, bs="tp", k=4), data=charice )
summary(mgcv.gl02)

plot.gam(mgcv.gl02, residuals=T, se=T, shade=T, pch=16, cex=0.5)
vis.gam(mgcv.gl02, plot.type="persp", type="response", color="gray", theta=-60)
gam.check(mgcv.gl02)
anova(mgcv.null.gl, mgcv.gl02,test="F")

# dc13 and ch4 interaction
mgcv.gl03 <- gam(globe ~ te(dc13*ch4, bs="tp", k=4) + te(ch4, bs="tp", k=4), data=charice )
summary(mgcv.gl03)

plot.gam(mgcv.gl03, residuals=T, se=T, shade=T, pch=16, cex=0.5)
vis.gam(mgcv.gl03, plot.type="persp", type="response", color="gray", theta=-60)
gam.check(mgcv.gl02)
anova(mgcv.null.gl, mgcv.gl02,test="F")
