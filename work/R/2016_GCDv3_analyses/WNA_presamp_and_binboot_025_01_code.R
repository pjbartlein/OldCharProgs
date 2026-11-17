datapath <- "E:/Projects/GPWG/work/pastCent/data/Rda"
dataname <- "WNA-zt2k.Rda"
basename <- "zt2k"
binbeg <- -60
binend <- 24000
binstep <- 10
binname <- "bw10" # prebinning label
abw <- 20 # window half-width (smoothing parameter)
bootname <- "abw020"
fitname <- paste("WNA_bin_influx_",basename,"_",binname,"_",bootname,sep="")
nreps <- 500
agebinpath <- "e:/Projects/GPWG/work/pastCent/data/csv/"
agebinname <- paste(fitname,".csv",sep="")
agebinfile <- paste(agebinpath, agebinname, sep="")
binage <- seq(binbeg, binend, by=binstep)
head(binage); tail(binage)
yrbeg <- 1010; yrend <- 2010; yrstep=abw
YearCE <- as.integer(seq(yrbeg,yrend,by=yrstep))
head(YearCE); tail(YearCE)
pltYearCE <- YearCE + (yrstep/2)
min_yearCE <- yrbeg-(yrstep/2); max_YearCE <- yrend+(yrstep/2)
xmin <- 1500; xmax <- 2010; ymin1 <- -5.0; ymax1 <- 5.0; ymin2 <- -2.0; ymax2 <- 2.0
xlim=c(xmin,xmax); ylim1=c(ymin1,ymax1); ylim2 <- c(ymin2,ymax2)
xlab="Year CE"; xminortick <- 10
load(paste(datapath,dataname,sep="/"))
names(all)
str(all)
nsites <- dim(all$Age)[2]
maxbins <- length(binage)
xp <- matrix(NA, ncol = nsites, nrow = maxbins)
yp <- matrix(NA, ncol = nsites, nrow = maxbins)
wp <- matrix(NA, ncol = nsites, nrow = maxbins)
# loop over sites, prebinning and stacking results
for (i in 1:nsites) {
  # check to see if all data values are NA, if so, skip this site
  sitedata <- as.data.frame(cbind(all$Age[,i],all$TransData[,i]))
  names(sitedata)=c("age","zt")
  sitedata <- na.omit(sitedata)
  nonmiss <- length(sitedata$age)
  
  # binning
  if (nonmiss > 0) {
    
    # add a column of 1's for counting
    sitedata$one <- rep(1,nonmiss)
    
    # this definition of bin number seems to match that implicit in presample.f90
    binnum <- as.integer(ceiling((sitedata$age - binbeg - (binstep/2))/binstep)) + 1
    #head(cbind(sitedata$age,sitedata$zt,binnum,binage[binnum]), 20)
    #tail(cbind(sitedata$age,sitedata$zt,binnum,binage[binnum]), 20)
    
    # get average zt values (and average ages) for the data in each bin
    binave <- tapply(sitedata$zt, binnum, mean)
    binaveage <- tapply(sitedata$age, binnum, mean)
    bincount <- tapply(sitedata$one, binnum, sum)
    #head(cbind(binave,binaveage,bincount), 20)
    #tail(cbind(binave,binaveage,bincount), 20)
    
  }
  
  # save data 
  nbins <- length(binave)
  xp[1:nbins,i] <- 1950.0 - binaveage
  yp[1:nbins,i] <- binave
  wp[1:nbins,i] <- rep(1,nbins)
  head(cbind(xp[,i],yp[,i],wp[,i]), 20)
  tail(cbind(xp[,i],yp[,i],wp[,i]), 20)
}
x <- as.vector(xp)
y <- as.vector(yp)
w <- as.vector(wp)
lfdata <- data.frame(x,y,w)
lfdata <- na.omit(lfdata)
lfdata <- lfdata[lfdata$x > min_yearCE & lfdata$x <= max_YearCE, ]
x <- lfdata$x; y <- lfdata$y; w <- lfdata$w
binnum <- as.integer(ceiling((x-yrbeg-(yrstep/2.0))/yrstep))+1
#temp <- cbind(x,y,binnum)
binave <- tapply(y, binnum, mean)
binaveage <- tapply(x, binnum, mean)
bin_fit_all <- binave

# plot the bin averages -- note pltYearCE offset to center the plot steps
plot(NULL, ylim=ylim1, xlim=xlim, ylab="Z-scores Tranformed Influx", xlab=xlab, 
  cex.sub=0.8, sub=fitname, type="n")
axis(side = 1, at = seq(xmin, xmax, by = xminortick), labels = FALSE, tcl = -.25) 
axis(side = 1, at = seq(xmin, xmax, by = xminortick*5), labels = FALSE, tcl = -.5) 
points(x,y, pch=16, cex=0.5, col=rgb(0.5,0.5,0.5,0.70))
lines(bin_fit_all ~ pltYearCE, type="S", col="red", lwd=2)
#points(bin_fit_all ~ YearCE, col="black", pch=16, lwd=2)
#temp2 <- cbind(YearCE, bin_fit_all)
# Step 1 -- Set up to plot individual replications
plot(NULL, NULL, ylim=ylim2, xlim=xlim, ylab="Z-scores Transformed Influx", xlab=xlab, 
  cex.sub=0.8, sub=fitname, type="n")
axis(side = 1, at = seq(xmin, xmax, by = xminortick), labels = FALSE, tcl = -.25)
axis(side = 1, at = seq(xmin, xmax, by = xminortick*5), labels = FALSE, tcl = -.5) 
yfit <- matrix(NA, nrow=length(YearCE), ncol=nreps)
set.seed(10) # do this to get the same sequence of random samples for each run
sites <- seq(1:nsites)

# Step 2 -- Do the bootstrap iterations, and plot each age-bin curve
ptm <- proc.time() # time the loop
for (i in 1:nreps) {
  #print(i)
  randsitenum <- sample(seq(1:nsites), nsites, replace=TRUE)
  # print(head(randsitenum))
  x <- as.vector(xp[,randsitenum])
  y <- as.vector(yp[,randsitenum])
  lfdata <- data.frame(x,y)
  lfdata <- na.omit(lfdata)
  lfdata <- lfdata[lfdata$x > min_yearCE & lfdata$x <= max_YearCE, ]
  x <- lfdata$x; y <- lfdata$y

  binnum <- as.integer(ceiling((x-yrbeg-(yrstep/2.0))/yrstep))+1
  binave <- tapply(y, binnum, mean)
  binaveage <- tapply(x, binnum, mean)
  bin_fit <- binave
  yfit[,i] <- bin_fit
  pltYearCE <- YearCE + (yrstep/2)
  lines(pltYearCE, yfit[,i], lwd=2, type="S", col=rgb(0.5,0.5,0.5,0.10))
}
proc.time() - ptm # how long?

# Step 3 -- Plot the unresampled (initial) area averages
lines(pltYearCE, bin_fit_all, type="S", lwd=2, col="red")

# Step 4 -- Find and add bootstrap CIs
yfit975 <- apply(yfit, 1, function(x) quantile(x,prob=0.975, na.rm=T))
yfit025 <- apply(yfit, 1, function(x) quantile(x,prob=0.025, na.rm=T))
yfit50 <- apply(yfit, 1, function(x) quantile(x,prob=0.500, na.rm=T))
lines(pltYearCE, yfit975, type="S", lwd=1, col="red")
lines(pltYearCE, yfit025, type="S", lwd=1, col="red")

outdata <- data.frame(YearCE,bin_fit_all,yfit025,yfit975,yfit50)
colnames(outdata) <- c("YearCE", fitname, "yfit025","yfit975","yfit50")
write.table (outdata, agebinfile, row.names=FALSE, sep=",")
