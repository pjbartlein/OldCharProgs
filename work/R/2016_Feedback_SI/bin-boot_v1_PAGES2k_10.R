# age-bin averages of temperature
# bootstrap-by-sample, stratfied by bins confidence intervals

# paths for input and output .csv files -- modify as appropriate
datapath <- "e:/Projects/GPWG/work/feedback_old/data/ice_core/"
outpath <- datapath

# number of bootstrap replications
nreps <- 200

# age bin centers
abw <-  20 # 100 #  50 #  
abinbeg <- -60 # -50 # -50 # 
abinend <- 2200

dataname <- "PAGES2k_tave"

# plot output 
plotout <- "pdf" # "screen" #  "png" # 

# plotting set up
xmin <- 0; xmax <- 2020; ymin1 <- -1.0; ymax1 <- 1.0; ymin2 <- -1.0; ymax2 <- 1.0
xlim=c(xmin,xmax); ylim1=c(ymin1,ymax1); ylim2 <- c(ymin2,ymax2)
xlab="Year CE"; xminortick <- 50
ylab="PAGES2k temperature"

# no changes below here

inputname <- paste(dataname, ".csv", sep="")
inputname

# curve (output) file
curvefile <- paste(dataname,"_binboot_v1_bw",as.character(abw),"_", as.character(nreps), ".csv", sep="")
print(curvefile)

# age bin centers
abinage <- seq(abinbeg, abinend, by=abw)
abinage
YearCE <- 1950 - abinage
YearCE

# generate YearCE with a 1/2 time-step (abw) offset so step-plot type "s" registers correctly
pltYearCE <- YearCE  +  (abw/2)
pltYearCE; length(pltYearCE)

# matrices for data and fitted values
maxreps <- 1000
yfit <- matrix(NA, nrow=length(abinage), ncol=maxreps)

# array for bootstrap results
min_age <- abinbeg-(abw/2); max_age <- abinend #+(abw/2)
nbins <- length(abinage)
yfit <- matrix(NA, nrow=nbins, ncol=maxreps)
 

# .pdf plot of bootstrap iterations
if (plotout == "pdf") {
  pdffile <- paste(dataname,"_binboot_v1_bw",as.character(abw),"_", as.character(nreps),".pdf", sep="")
  print(pdffile)
}
# .png plot of bootstrap iterations
if (plotout == "png") {
  pngfile <- paste(dataname,"_binboot_v1_bw",as.character(abw),"_", as.character(nreps),".png", sep="")
  print(pngfile)
}


# read data
inputfile <- paste(datapath, inputname, sep="")
print(inputfile)
indata <- read.csv(inputfile)
head(indata); tail(indata, 20)

age <- indata$Age_yrBP1950 # Age_yrBP1950, tave
depvar <- indata$tave 

# trim samples to age range
depvar <- depvar[age < abinend+abw] 
age <- age[age < abinend+abw] 
nrecs <- length(age)


# reshape matrices into vectors 
x <- as.vector(age)
y <- as.vector(depvar)
lfdata <- data.frame(x,y)
lfdata <- na.omit(lfdata)
x <- lfdata$x; y <- lfdata$y
one <- rep(1,length(x))

# average value for each age bin

binnum <- as.integer(ceiling((x-abinbeg-(abw/2.0))/abw))+1
binave <- tapply(y, binnum, mean)
binaveage <- tapply(x, binnum, mean)
bincount <- tapply(one, binnum, sum)
binsubs <- as.integer(unlist(dimnames(binave)))
bin_fit_all <- rep(NA, length(abinage))
bin_fit_all[binsubs] <- binave
ninbin <- rep(NA, length(abinage))
ninbin[binsubs] <- bincount
                 
#print(binnum); print(length(binnum))
#print(bincount); length(bincount)
#print(binave); print(length(binave))
#print(binaveage); print(length(binaveage))
#print(binsubs); print(length(binsubs))
#print(bin_fit_all); print(length(bin_fit_all))
#print(abinage); print(length(abinage)) 
#print(ninbin); print(length(ninbin))
#print(pltYearCE); print(length(pltYearCE)) 
#print(cbind(abinage[binsubs], bincount, binave, binsubs, 1950-abinage[binsubs]))

#print(cbind(abinage, ninbin, bin_fit_all, 1950-abinage))

#cbind(abinage, bin_fit_all, ninbin)
#head(cbind(1950 - binaveage, pltYearCE[binsubs], abinage[binsubs], bin_fit_all[binsubs]))
#tail(cbind(1950 - binaveage, pltYearCE[binsubs], abinage[binsubs], bin_fit_all[binsubs]))

# Bootstrap samples

# Step 1 -- Set up to plot individual replications
if (plotout == "pdf") {pdf(file=paste(outpath,pdffile,sep=""))}
if (plotout == "png") {png(file=paste(outpath,pngfile,sep=""), res=150)}
plot(NULL, ylim=ylim1, xlim=xlim, ylab=ylab, xlab=xlab, cex.sub=0.8, sub=curvefile, type="n")
axis(side = 1, at = seq(xmin, xmax, by = xminortick*5), labels = FALSE, tck=1.0, fg="gray70")
axis(side = 1, at = seq(xmin-xminortick, xmax+xminortick, by = xminortick), labels = FALSE, tcl = -.25)
axis(side = 1, at = seq(xmin, xmax, by = xminortick*5), labels = FALSE, tcl = -.5) 
 

points(1950-indata$Age_yrBP1950, indata$tave, pch=16, cex=0.3, col="lightblue")

# for debugging step plots
# plot the bin averages -- note pltYearCE offset to center the plot steps
# points(1950 - x, y, pch=16, cex=0.5, col=rgb(0.5,0.5,0.5,0.70))
# lines(bin_fit_all ~ pltYearCE, type="s", col="red", lwd=2)
# points(1950 - abinage, bin_fit_all, col="blue", pch=16, cex=0.5)

set.seed(10) # do this to get the same sequence of random samples for each run

# Step 2 -- Do the bootstrap iterations, and plot each age-bin curve
ptm <- proc.time() # time the loop
for (i in 1:nreps) {
  print(i)
  randrecnum <- sample(seq(1:nrecs), nrecs, replace=TRUE)
  # print(head(randrecnum))
  x <- as.vector(age[randrecnum])
  y <- as.vector(depvar[randrecnum])
  lfdata <- data.frame(x,y)
  lfdata <- na.omit(lfdata)
  lfdata <- lfdata[lfdata$x >= min_age & lfdata$x < max_age, ]
  lfdata_boot <- lfdata[randrecnum, ]
  
  #x_bs <- NA; y_bs <- NA
  #lfdata_boot <- data.frame(x_bs, y_bs)
  #names(lfdata_boot) <- c("x", "y")
  #for (j in 1:nbins) {
  #  #j <- 1
  #  if (!is.na(ninbin[j])) {
  #    lfdata_bs <- lfdata[binnum == j, ]
  #    randrecnum <- sample(seq(1:ninbin[j]), replace=TRUE)
  #    lfdata_boot <- rbind(lfdata_boot, lfdata_bs[randrecnum,])
  #  }
  #}
  
  lfdata_boot <- na.omit(lfdata_boot)
  x <- lfdata_boot$x; y <- lfdata_boot$y
  one <- rep(1,length(x))
  
  binnum_bs <- rep(NA, nrecs)
  binnum_bs <- as.integer(ceiling((x-abinbeg-(abw/2.0))/abw))+1
  binave_bs <- rep(NA,nbins)
  binave_bs <- tapply(y, binnum_bs, mean)
  binaveage_bs <- rep(NA,nbins)
  binaveage_bs <- tapply(x, binnum_bs, mean)
  bincount_bs <- tapply(one, binnum_bs, sum)
  binsubs_bs <- as.integer(unlist(dimnames(binave_bs)))
  bin_fit_all_bs <- rep(NA,nbins)
  bin_fit_all_bs[binsubs_bs] <- binave_bs
  ninbin <- rep(NA, length(abinage))
  ninbin[binsubs_bs] <- bincount_bs
  yfit[,i] <- bin_fit_all_bs
  
  #print(binnum_bs); print(length(binnum_bs))
  #print(abinage); length(abinage)
  #print(binave_bs); print(length(binave_bs))
  #print(binsubs_bs); print(length(binsubs_bs))
  
  #print(cbind(abinage[binsubs_bs], bincount_bs, binave_bs, binsubs_bs, 1950-abinage[binsubs_bs]))
  #print(cbind(abinage, ninbin, bin_fit_all_bs, yfit[,i], 1950-abinage))
  #print(tail(cbind(1950-sort(x), y[order(x)], binnum_bs[order(x)],1950- abinage[binnum_bs[order(x)]]), 30))
  
  # print(c(x[binnum_bs == 39], y[binnum_bs == 39]))
  # print(c(bin_fit_all_bs[39], i, yfit[39,i]))
                      
  segments(pltYearCE-1.9*(abw/2), yfit[,i], pltYearCE-0.1*(abw/2), yfit[,i], lwd=0.6, col=rgb(0.3,0.3,0.3,0.25))
  #points(pltYearCE, yfit[,1])
}
proc.time() - ptm # how long?
#tail(cbind(pltYearCE, yfit[,1], yfit[,2], yfit[,3], yfit[,4], yfit[,i-2], yfit[,i-1],yfit[,i], ninbin), 20)

# Step 3 -- Plot the unresampled (initial) area averages
lines(pltYearCE[binsubs], bin_fit_all[binsubs], type="s", lwd=1, col="red")
lastbin <- binsubs[length(binsubs)]
segments(pltYearCE[lastbin]-abw, bin_fit_all[lastbin], pltYearCE[lastbin], bin_fit_all[lastbin], lwd=1, col="red")
#points(pltYearCE[binsubs]-(abw/2), bin_fit_all[binsubs], pch=16, cex=0.8, col="blue")

# Step 4 -- Find and add bootstrap CIs
yfit975 <- apply(yfit, 1, function(x) quantile(x,prob=0.975, na.rm=T))
yfit025 <- apply(yfit, 1, function(x) quantile(x,prob=0.025, na.rm=T))
yfit50 <- apply(yfit, 1, function(x) quantile(x,prob=0.500, na.rm=T))

#lines(pltYearCE, yfit975, type="s", lwd=0.5, col="red")
segments(pltYearCE-1.9*(abw/2), yfit975, pltYearCE-0.1*(abw/2), yfit975, lwd=0.6, col="red")
segments(pltYearCE[nbins]-abw, yfit975[nbins], pltYearCE[nbins], yfit975[nbins], lwd=1, col="red")

#lines(pltYearCE, yfit025, type="s", lwd=0.5, col="red")
segments(pltYearCE-1.9*(abw/2), yfit025, pltYearCE-0.1*(abw/2), yfit025, lwd=0.6, col="red")
segments(pltYearCE[nbins]-abw, yfit025[nbins], pltYearCE[nbins], yfit025[nbins], lwd=1, col="red")

if (plotout == "pdf") {dev.off()}
if (plotout == "png") {dev.off()}
curveout <- data.frame(cbind(abinage, 1950-abinage, bin_fit_all, yfit975, yfit025, ninbin))
colnames(curveout) <- c("age", "YearCE", "bin_ave", "cu95", "cl95", "ninbin")
outputfile <- paste(outpath, curvefile, sep="")
write.table(curveout, outputfile, col.names=TRUE, row.names=FALSE, sep=",")




