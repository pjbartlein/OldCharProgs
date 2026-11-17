# age-bin averages of Mann et al. temperature
# bootstrap-by-sample

# paths for input and output .csv files -- modify as appropriate
datapath <- "e:/Projects/GPWG/work/feedback/data/ice_core/"
outpath <- datapath

library(locfit)

# locfit (half) window-width parameter
# note:  hw should be an integer multiple of pbw
hw <- 50 # bandwidth (smoothing parameter)
hw_char <- as.character(hw)
if (hw < 100) hw_char <- paste("0", hw_char, sep="")

# number of bootstrap samples/replications
nreps <- 200

# target ages for fitted values
targstep <- 50
targbeg <- -50
targend <- 1450

dataname <- "aaves_tmp_land"

# plot output 
plotout <- "pdf" # "screen" #  "png" # 

# plotting set up
xmin <- 0; xmax <- 2020; ymin1 <- -1.0; ymax1 <- 1.0; ymin2 <- -1.0; ymax2 <- 1.0
xlim=c(xmin,xmax); ylim1=c(ymin1,ymax1); ylim2 <- c(ymin2,ymax2)
xlab="Year CE"; xminortick <- 50
ylab="Mann et al. temperature"

# no changes below here

# target values
targage <- seq(targbeg,targend,targstep)
targage
YearCE <- 1950 - targbeg
YearCE

inputname <- paste(dataname, ".csv", sep="")
inputname

# curve (output) file
curvename <- paste(dataname,"_lwrboot_hw",hw_char,"_", as.character(nreps), sep="")
curvefile <- paste(curvename, ".csv", sep="")
print(curvefile)

# generate YearCE with a 1/2 time-step (abw) offset so step-plot type "s" registers correctly
pltYearCE <- YearCE  +  (targstep/2)
pltYearCE; length(pltYearCE)

# matrices for data and fitted values
maxreps <- 1000
min_age <- targbeg-(targstep/2); max_age <- targend + (targstep/2)
yfit <- matrix(NA, nrow=length(targage), ncol=maxreps)

# array for bootstrap results
targage <- seq(targbeg,targend,targstep)
targage.df <- data.frame(x=targage)
lowage <- targage - hw; highage <- targage + hw
ntarg <- length(targage)
yfit <- matrix(NA, nrow=length(targage.df$x), ncol=maxreps)
 

# .pdf plot of bootstrap iterations
if (plotout == "pdf") {
  pdffile <- paste(dataname,"_lwrboot_hw",hw_char,"_", as.character(nreps),".pdf", sep="")
  print(pdffile)
}
# .png plot of bootstrap iterations
if (plotout == "png") {
  pngfile <- paste(ddataname,"_lwrboot_hw",hw_char,"_", as.character(nreps),".png", sep="")
  print(pngfile)
}


# read data
inputfile <- paste(datapath, inputname, sep="")
print(inputfile)
indata <- read.csv(inputfile)
head(indata); tail(indata, 20)

age <- indata$Age_yrBP1950 # Age_yrBP1950
depvar <- indata$GLOBE

# trim samples to age range
depvar <- depvar[age < targend+targstep] 
age <- age[age < targend+targstep] 
nrecs <- length(age)


# reshape matrices into vectors 
x <- as.vector(age)
y <- as.vector(depvar)
lfdata <- data.frame(x,y)
lfdata <- na.omit(lfdata)
x <- lfdata$x; y <- lfdata$y
one <- rep(1,length(x))

# 2. locfit
# initial fit, unresampled (i.e. all) data
ptm <- proc.time()
loc01 <- locfit(y ~ lp(x, deg=1, h=hw), maxk=800, maxit=20, family="qrgauss")
summary(loc01)

# 3. get  fitted values
pred01 <- predict(loc01, newdata=targage.df, se.fit=TRUE)
loc01_fit <- data.frame(targage.df$x, pred01$fit)
fitname <- paste("locfit_",as.character(hw), sep="")
colnames(loc01_fit) <- c("age", fitname)
head(loc01_fit)
proc.time() - ptm
ptm <- proc.time()

# Bootstrap samples

# Step 1 -- Set up to plot individual replications
if (plotout == "pdf") {pdf(file=paste(outpath,pdffile,sep=""))}
if (plotout == "png") {png(file=paste(outpath,pngfile,sep=""), res=150)}
plot(NULL, NULL, ylim=ylim2, xlim=xlim, ylab=ylab, xlab=xlab, cex.sub=0.8, sub=curvename, type="n")
axis(side = 1, at = seq(xmin, xmax, by = xminortick*5), labels = FALSE, tck=1.0, fg="gray70")
axis(side = 1, at = seq(xmin-xminortick, xmax+xminortick, by = xminortick), labels = FALSE, tcl = -.25)
axis(side = 1, at = seq(xmin, xmax, by = xminortick*5), labels = FALSE, tcl = -.5) 

points(1950-indata$Age_yrBP1950, indata$GLOBE, pch=16, cex=0.3, col="lightblue")

# Step 2 -- Do the bootstrap iterations, and plot each composite curve
set.seed(10) # do this to get the same sequence of random samples for each run

nerr <- 0
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

  x <- lfdata_boot$x; y <- lfdata_boot$y
  locboot <- locfit(y ~ lp(x, deg=1, h=hw), maxk=400, maxit=20, debug=0, family="qrgauss")
  result <- try(predict(locboot, newdata=targage.df, se.fit=TRUE), silent=TRUE)
  if (class(result) != "try-error") {
    predboot <- predict(locboot, newdata=targage.df, se.fit=TRUE)
    yfit[,i] <- predboot$fit
    lines(1950 - targage.df$x, yfit[,i], lwd=2, col=rgb(0.5,0.5,0.5,0.10))
  } else {
    print ("NA's produced")
    i <- i-1
    nerr <- nerr + 1
  }
  if (i %% 10 == 0) {print(i)}
}
print(nerr)

# Step 3 -- Plot the unresampled (initial) fit
fitname <- paste("locfit_",as.character(hw), sep="")
colnames(loc01_fit) <- c("age", fitname)
lines(1950 - loc01_fit[,1], loc01_fit[,2], lwd=2, col="red")
#points(1950 - loc01_fit[,1], loc01_fit[,2], pch=16, cex=0.5, col="blue")

# Step 4 -- Find and add bootstrap CIs
yfit975 <- apply(yfit, 1, function(x) quantile(x,prob=0.975, na.rm=T))
yfit025 <- apply(yfit, 1, function(x) quantile(x,prob=0.025, na.rm=T))
lines(1950 - targage.df$x, yfit975, lwd=1, col="red")
lines(1950 - targage.df$x, yfit025, lwd=1, col="red")

if (plotout == "pdf") {dev.off()}
if (plotout == "png") {dev.off()}
curveout <- data.frame(cbind(targage, 1950-targage, loc01_fit[,2], yfit975, yfit025))
colnames(curveout) <- c("age", "YearCE", "locfit", "cu95", "cl95")
outputfile <- paste(outpath, curvefile, sep="")
write.table(curveout, outputfile, col.names=TRUE, row.names=FALSE, sep=",")




