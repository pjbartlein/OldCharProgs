setwd("/Users/Olivier/Documents/paleofiRe_Paper/")

library(paleofiRe)

## Sites selection process
ID=pfSiteSel(Latlim=c(44,90), Longlim=c(-90,-50), DateInt=2500)

sumID=summary(ID)
sites_inc=sumID$ID_SITE[sumID$NUM_SAMP>=20]
ID=pfSiteSel(ID=sites_inc)

## Figure 1
postscript("Fig1_n.eps",horizontal=FALSE,width=10, height=12, paper = "special")
par(mfrow=c(2,1))
plot(ID,zoom="world")
plot(ID, type="Chronology")
dev.off()

## Data transformation
TR1=pfTransform(ID, method=c("MinMax", "Box-Cox", "ZScore"))
TR2=pfTransform(ID, method=c("Box-Cox","RunMean","ZScore"), type="JohnDraper", 
               RunWidth=500,Interpolate=TRUE, Age=seq(-50,11500,by=50))

## Data compositing
COMP1=pfComposite(TR1, binning=TRUE, bins=seq(0,11000,250))
COMP2=pfComposite(TR2, binning=FALSE, nboot=1000)

summary(COMP2)

## Figure 2

postscript("Fig2_n.eps",horizontal=FALSE,width=6, height=8, paper = "special")
plot(COMP2,conf=c(0.005,0.025,0.975,0.995))
dev.off()

## Circular block bootstrapp
sea1=SEA(COMP2,conf=c(0.005,0.025,0.975,0.995),nboot=100)

# Figure 3
plot(sea1)



### TEST
## pfKrige usage not presented in paper
## Test: compare 0 6000 bp...

bp6=pfKrige(TR1,AgeLim=c(6500,5500),gridkm=50)
p2=plot(bp6)
bp0=pfKrige(TR1,AgeLim=c(0,500),gridkm=50)
p1=plot(bp0)


require(gridExtra)
grid.arrange(p1+ ggtitle("0-500 BP"), p2+ ggtitle("5500-6500 BP"), nrow=2)

## 


