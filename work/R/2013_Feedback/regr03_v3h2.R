library(visreg)


## ------------------------------------------------------------------------
# read the charcoal and ice-core data
datapath <- "e:/Projects/GPWG/work/feedback/data/ice_core/"
datafile <- "charice_binned_01.csv"
input <- read.csv(paste(datapath,datafile,sep=""))
input$year_CE <- 1950-input$targ_age
summary(input)


## ------------------------------------------------------------------------
# function for plotting two series at a time
plotseries <- function(title, dataset, xvarin, yvar1in, yvar2in) {
  xvar <- dataset[,xvarin]
  yvar1 <- dataset[,yvar1in]; yvar2 <- dataset[,yvar2in]
  opar <- par(mar=c(5,4,4,5)+0.1) # space for second axis
  plot(yvar1 ~ xvar, pch=16, xlab=xvarin, ylab=yvar1in, col="red", data=dataset,
    main=title)
  par(new=T)          # second plot is going to get added to first
  plot(yvar2 ~ xvar, xlab="", pch=16, axes=F, ylab="", col="blue", data=dataset)
  axis(side=4); mtext(side=4,line=3.8,yvar2in)
  legend("bottom", c(yvar1in,yvar2in), pch=16, col=c("red","blue"))
  par(opar)  # restore plot parameters
}


## plot all data
charice <- input

## ----fig.width=8, fig.height=6-------------------------------------------
plotseries("global charcoal and (adjusted) northern & southern d13C", charice, "year_CE", "globe", "nsd13c_adj")


## ------------------------------------------------------------------------
# select the data between 0 and 1700 CE
charice <- input[input$targ_age >= 250.0 & input$targ_age <= 1950,]
# number of observations after trimming
dim(charice)[1]  # number of observations
#summary(charice)

## ----fig.width=8, fig.height=6-------------------------------------------
plotseries("global charcoal and (adjusted) northern & southern d13C", charice, "year_CE", "globe", "nsd13c_adj")


## ------------------------------------------------------------------------
plot (nsd13c_adj ~ globe, data=charice, pch=16)
text (nsd13c_adj ~ globe, labels=as.character(year_CE), data=charice, pos=3, cex=.75)



## ------------------------------------------------------------------------
charice2 <- input[input$targ_age >= 350.0 & input$targ_age <= 1950,]
dim(charice2)[1]  # number of observations
charice3 <- input[input$targ_age >= 350.0 & input$targ_age <= 1850,]
dim(charice3)[1]  # number of observations


## ----fig.width=8, fig.height=6-------------------------------------------
plotseries("global charcoal and (adjusted) northern & southern d13C", charice2, "year_CE", "globe", "nsd13c_adj")

## ------------------------------------------------------------------------
plot (nsd13c_adj ~ globe, data=charice2, pch=16)
text (nsd13c_adj ~ globe, labels=as.character(year_CE), data=charice2, pos=3, cex=.75)


## ----fig.width=8, fig.height=6-------------------------------------------
plotseries("global charcoal and (adjusted) northern & southern d13C", charice3, "year_CE", "globe", "nsd13c_adj")

## ------------------------------------------------------------------------
plot (nsd13c_adj ~ globe, data=charice2, pch=16)
text (nsd13c_adj ~ globe, labels=as.character(year_CE), data=charice3, pos=3, cex=.75)



## ------------------------------------------------------------------------
# regression -- d13C as a function of charcoal
lm.nsdc13_adj.globe <- lm(nsd13c_adj ~ globe, data=charice)
summary(lm.nsdc13_adj.globe)
# some diagnostic plots
hist(lm.nsdc13_adj.globe$residuals, breaks=20)
qqnorm(lm.nsdc13_adj.globe$residuals)
qqline(lm.nsdc13_adj.globe$residuals)
plot(lm.nsdc13_adj.globe$residuals ~ charice$year_CE)


## charcoal as a function of d13C
plot (globe ~nsd13c_adj, data=charice, pch=16)
text (globe ~nsd13c_adj, labels=as.character(year_CE), data=charice, pos=3, cex=.75)


## ------------------------------------------------------------------------
# "calibration" regression with CH4 and dC13 * CH4 product term as predictors
lm.globe.nsd13c_adj.calib <- lm(globe ~ I(nsd13c_adj*CH4) + CH4, data=charice)
summary(lm.globe.nsd13c_adj.calib)


## ------------------------------------------------------------------------
hist(lm.globe.nsd13c_adj.calib$residuals, breaks=20)
qqnorm(lm.globe.nsd13c_adj.calib$residuals)
qqline(lm.globe.nsd13c_adj.calib$residuals)
plot(lm.globe.nsd13c_adj.calib$residuals ~ charice$targ_age, xlim=c(2000,0))


##
visreg(lm.globe.nsd13c_adj.calib, points=list(pch=16, cex=1))

visreg2d(lm.globe.nsd13c_adj.calib, xvar="nsd13c_adj", yvar= "CH4")