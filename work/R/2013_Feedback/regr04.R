# regression of influx on dc13 and ch4

# function for plotting two series at a time
plotseries <- function(title, dataset, xvarin, yvar1in, yvar2in) {
  print(title)
  xvar <- dataset[,xvarin]
  yvar1 <- dataset[,yvar1in]; yvar2 <- dataset[,yvar2in]
  opar <- par(mar=c(5,4,4,5)+0.1) # space for second axis
  plot(yvar1 ~ xvar, pch=16, xlab=xvarin, ylab=yvar1in, col="red", data=dataset,
    main=title)
  par(new=T)          # second plot is going to get added to first
  plot(yvar2 ~ xvar, xlab="", pch=16, axes=F, ylab="", col="blue", data=dataset)
  axis(side=4); mtext(side=4,line=3.8,yvar2in)
  legend("bottom", c(yvar1in,yvar2in), pch=16, col=c("red","blue"))
  par(opar)  # restore plot par
}


# read the charcoal and ice-core data
datapath <- "e:/Projects/GPWG/work/feedback/data/ice_core/"
datafile <- "charice_02.csv"
input <- read.csv(paste(datapath,datafile,sep=""))
input$year_CE <- 1950-input$targ_age

# select the data between 0 and 1700 CE
charice <- input[input$targ_age >= 250.0 & input$targ_age <= 1950,]
summary(charice)
dim(charice)

# plot a charcoal and dc13 curve
title <- "globe and northern d13C"
plotseries(title, charice, "year", "globe", "nd13c")

# scatter plot
plot (nd13c ~ globe, data=charice, pch=16)
text (nd13c ~ globe, labels=as.character(year_CE), data=charice, pos=3, cex=.5)

# regression -- d13C as a function of charcoal
lm.ndc13.globe <- lm(nd13c ~ globe, data=charice)
summary(lm.ndc13.globe)
hist(lm.ndc13.globe$residuals, breaks=20)
plot(lm.ndc13.globe$residuals ~ charice$year)

# "calibration" regression with dc13 and ch4 product term
lm.globe.nd13c.calib <- lm(globe ~ I(nd13c*ch4) + ch4, data=charice)
summary(lm.globe.nd13c.calib)
hist(lm.globe.nd13c.calib$residuals, breaks=20)
plot(lm.globe.nd13c.calib$residuals ~ charice$targ_age, xlim=c(2000,0))

plotseries(title, charice, "year", "nhemis", "nsd13c_adj")

# scatter plot
plot (nsd13c_adj ~ nhemis, data=charice, pch=16)
text (nsd13c_adj ~ nhemis, labels=as.character(year), data=charice, pos=3, cex=.5)
