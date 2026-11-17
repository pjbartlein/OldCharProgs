# regression of influx on dc13 and ch4

datapath <- "e:/Projects/GPWG/work/feedback/data/ice_core/"
datafile <- "charice_100.csv"

input <- read.csv(paste(datapath,datafile,sep=""))
charice <- input[input$targ_age >= 250.0 & input$targ_age <= 2000,]
summary(charice)

# plot a charcoal and dc13 curve
opar <- par(mar=c(5,4,4,5)+0.1) # space for second axis
plot(nextr ~ targ_age, pch=16, col="red", data=charice)  # first plot
par(new=T)          # second plot is going to get added to first
plot(nd13c ~ targ_age, pch=16, axes=F, ylab="", col="gray", data=charice)  # second plot
axis(side=4)   # add axis
mtext(side=4,line=3.8,"dc13")   # add label
par(opar)  # restore plot par


# regression with dc13 and ch4 product term
lm.globe <- lm(globe ~ I(dc13*ch4) + ch4, data=charice)
summary(lm.globe)
plot(lm.globe, which=c(1,2,4))
hist(lm.globe$residuals, breaks=20)

