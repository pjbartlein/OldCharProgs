# composite curve via the locfit package
# bootstrap-by-site confidence intervals

# names
queryname <- "v3i" # 
basename <- "nt2kb" # "zt2kb" # 

# paths for input and output .csv files -- modify as appropriate
datapath <- "/Projects/GPWG/GPWGv3/GCDv3Data/v3i/"
sitelistpath <- "/Projects/GPWG/GPWGv3/GCDv3Data/v3i/v3i_sitelists/"
sitelist <- "v3i_nsa_globe"
outpath <- "/Projects/GPWG/GPWGv3/GCDv3Data/v3i/v3i_curves/"

# prebinning bin width
pbw <- 10

pbw_char <- as.character(pbw)
if (pbw < 10) pbw_char <- paste("0", pbw_char, sep="")

# prebinning file name
binname <- paste("bw", pbw_char, sep="")

# presampled/binned files
csvpath <- "/Projects/GPWG/GPWGv3/GCDv3Data/v3i/v3i_presamp_csv/"
csvname <- paste("_presamp_influx_",basename,"_",binname,".csv", sep="")

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
targend <- 2200

# array sizes
maxrecs <- 2000
maxreps <- 1000

# plotting set up
xmin <- 0; xmax <- 2020; ymin1 <- -1.0; ymax1 <- 1.0; ymin2 <- -1.0; ymax2 <- 1.0
xlim=c(xmin,xmax); ylim1=c(ymin1,ymax1); ylim2 <- c(ymin2,ymax2)
xlab="Year CE"; xminortick <- 50
#ylab="Z-Scores of Transformed Influx"
ylab="Normalized Anomalies of Transformed Influx"

# plot output 
plotout <- "pdf" # "screen" # 

# no changes below here

# prebinning file name
binname <- paste("bw", pbw_char, sep="")

# site list file
sitelistfile <- paste(sitelistpath, sitelist, ".csv", sep="")
sitelistfile

# curve (output) path and file
curvecsvpath <- paste(datapath,queryname,"_curves/",sep="")

# if output folder does not exist, create it
dir.create(file.path(datapath, paste(queryname,"_curves/",sep="")), showWarnings=FALSE)
curvename <- paste(sitelist,"_locfit_",basename,"_",binname,"_hw",hw_char,"_",
  as.character(nreps), sep="")
curvefile <- paste(curvename, ".csv", sep="")
print(curvecsvpath)
print(curvename)
print(curvefile)

# .pdf plot of bootstrap iterations
if (plotout == "pdf") {
pdffile <- paste(curvename, ".pdf", sep="")
print(pdffile)
}

# read the list of sites
sites <- read.csv(sitelistfile)
head(sites)
ns <- length(sites[,1]) #length(sites$ID_SITE)
ns
# arrays for data and fitted values
age <- matrix(NA, ncol=ns, nrow=maxrecs)
influx <- matrix(NA, ncol=ns, nrow=maxrecs)
nsamples <- rep(0, maxrecs)
targage <- seq(targbeg,targend,targstep)
targage.df <- data.frame(x=targage)
lowage <- targage - hw; highage <- targage + hw
ntarg <- length(targage)
yfit <- matrix(NA, nrow=length(targage.df$x), ncol=maxreps)

# arrays for sample number and effective window span tracking
ndec <- matrix(0, ncol=ntarg, nrow=ns)
ndec_tot <- rep(0, ntarg)
xspan <- rep(0, ntarg)
ninwin <- matrix(0, ncol=ntarg, nrow=ns)
ninwin_tot <- rep(0, ntarg)

# read and store the presample (binned) files as matrices of ages and influx values
ii <- 0
for (i in 1:ns) {
  #i <- 1
  sitenum <- sites[i,1] # sites$ID_SITE[i]
  print(sitenum)
  siteidchar <- as.character(sitenum)
  if (sitenum >= 1) siteid <- paste("000", siteidchar, sep="")
  if (sitenum >= 10) siteid <- paste("00", siteidchar, sep="")
  if (sitenum >= 100) siteid <- paste("0", siteidchar, sep="")
  if (sitenum >= 1000) siteid <- paste(    siteidchar, sep="")
  
  inputfile <- paste(csvpath, siteid, csvname, sep="")
  print(inputfile)
  
  if (file.exists(inputfile)) {
    indata <- read.csv(inputfile)
    nsamp <-  length(indata$age) # 
    if (nsamp > 0) {
        ii <- ii+1
        age[1:nsamp,ii] <- indata$age # 
        influx[1:nsamp,ii] <- indata$norman # indata$zt # 
        nsamples[ii] <- nsamp
    }
  }
}
nsites <- ii

# number of sites with data
nsites

# trim samples to age range
influx[age >= targend+hw] <- NA
age[age >= targend+hw] <- NA

# censor abs(influx) values > 10
influx[abs(influx) >= 10] <- NA
age[abs(influx) >= 10] <- NA

# count number of sites that contributed to each fitted value
ptm <- proc.time()
for (i in 1:ntarg) {
  agemax <- -1e32; agemin <- 1e32
  for (j in 1:nsites) {
    for (k in 1:nsamples[j]) {
      if (!is.na(age[k,j])) {
        ii <- as.integer(ceiling((age[k,j]-targbeg)/targstep))+1
        #print (c(i,j,k,ii))
        if (ii > 0 && ii <= ntarg) {ndec[j,ii] = 1}
        if (age[k,j] >= targage[i]-hw && age[k,j] <= targage[i]+hw) {
          ninwin[j,i] = 1
          if (agemax <= age[k,j]) {agemax <- age[k,j]}
          if (agemin >= age[k,j]) {agemin <- age[k,j]}
        }
      }
    }
  }
  ndec_tot[i] <- sum(ndec[,i])
  ninwin_tot[i] <- sum(ninwin[,i])
  xspan[i] <- agemax - agemin
}
proc.time() - ptm
head(cbind(targage,ndec_tot,xspan,ninwin_tot))
tail(cbind(targage,ndec_tot,xspan,ninwin_tot))

ptm <- proc.time()
# 1. reshape matrices into vectors 
x <- as.vector(age)
y <- as.vector(influx)
lfdata <- data.frame(x,y)
lfdata <- na.omit(lfdata)
x <- lfdata$x; y <- lfdata$y

# 2. locfit
# initial fit, unresampled (i.e. all) data
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
if (plotout == "pdf") {pdf(file=paste(curvecsvpath,pdffile,sep=""))}
plot(NULL, NULL, ylim=ylim2, xlim=xlim, ylab=ylab, xlab=xlab, cex.sub=0.8, sub=curvename, type="n")
axis(side = 1, at = seq(xmin-xminortick, xmax+xminortick, by = xminortick), labels = FALSE, tcl = -.25)
axis(side = 1, at = seq(xmin, xmax, by = xminortick*5), labels = FALSE, tcl = -.40)

# Step 2 -- Do the bootstrap iterations, and plot each composite curve
set.seed(10) # do this to get the same sequence of random samples for each run

nerr <- 0
for (i in 1:nreps) {
  print(i)
  randsitenum <- sample(seq(1:nsites), nsites, replace=TRUE)
  # print(head(randsitenum))
  x <- as.vector(age[,randsitenum])
  y <- as.vector(influx[,randsitenum])
  lfdata <- data.frame(x,y)
  lfdata <- na.omit(lfdata)
  x <- lfdata$x; y <- lfdata$y
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
yfit95 <- apply(yfit, 1, function(x) quantile(x,prob=0.975, na.rm=T))
yfit05 <- apply(yfit, 1, function(x) quantile(x,prob=0.025, na.rm=T))
lines(1950 - targage.df$x, yfit95, lwd=1, col="red")
lines(1950 - targage.df$x, yfit05, lwd=1, col="red")

if (plotout == "pdf") {dev.off()}
curveout <- data.frame(cbind(targage.df$x, pred01$fit, yfit95, yfit05, ndec_tot, xspan, ninwin_tot))
colnames(curveout) <- c("age", "locfit", "cu95", "cl95", "nsites", "window", "ninwin")
outputfile <- paste(curvecsvpath, curvefile, sep="")
write.table(curveout, outputfile, col.names=TRUE, row.names=FALSE, sep=",")
proc.time() - ptm
warnings()
