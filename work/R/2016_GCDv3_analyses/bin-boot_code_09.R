# age-bin averages of influx data
# bootstrap-by-site confidence intervals

# names
queryname <- "v3plus" # "v3i" # 
basename <- "zt-lme"

# paths for input and output .csv files -- modify as appropriate
datapath <- "/Projects/GPWG/GPWGv3/GCDv3Data/v3plus/"
sitelistpath <- "/Projects/GPWG/GPWGv3/GCDv3Data/v3plus/v3plus_sitelists/"
sitelist <-  "v3i_wus" # "wus2_sites" # "v3plus_WUS" # 
outpath <- "/Projects/GPWG/GPWGv3/GCDv3Data/v3plus/v3plus_curves/"

# prebinning bin width
pbw <- 10

pbw_char <- as.character(pbw)
if (pbw < 10) pbw_char <- paste("0", pbw_char, sep="")

# number of bootstrap samples/replications
nreps <- 200

# age bin centers
abw <- 20
abinbeg <- -60 
abinend <- 1100

# array sizes
maxrecs <- 2000
maxreps <- 1000

# plotting set up
xmin <- 800; xmax <- 2020; ymin1 <- -5.0; ymax1 <- 5.0; ymin2 <- -2.0; ymax2 <- 2.0
xlim=c(xmin,xmax); ylim1=c(ymin1,ymax1); ylim2 <- c(ymin2,ymax2)
xlab="Year CE"; xminortick <- 20

# plot output 
plotout <- "pdf" # "screen" # 

# no changes below here

# prebinning file name
binname <- paste("bw", pbw_char, sep="")

# presampled/binned files
csvpath <- "/Projects/GPWG/GPWGv3/GCDv3Data/v3plus/v3plus_presamp_csv/"
csvname <- paste("_presamp_influx_",basename,"_",binname,".csv", sep="")

# site list file
sitelistfile <- paste(sitelistpath, sitelist, ".csv", sep="")
sitelistfile

# curve (output) path and file
curvecsvpath <- paste(datapath,queryname,"_curves/",sep="")

# if output folder does not exist, create it
dir.create(file.path(datapath, paste(queryname,"_curves/",sep="")), showWarnings=FALSE)
curvename <- paste(sitelist,"_binboot_",basename,"_",binname,"_abw",as.character(abw),"_",
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

# age bin centers
abinage <- seq(abinbeg, abinend, by=abw)
abinage
YearCE <- 1950 - abinage
YearCE

# generate YearCE with a 1/2 time-step (abw) offset so step-plot type "s" registers correctly
pltYearCE <- YearCE + (abw/2)
pltYearCE

#targage.df <- data.frame(x=targage)
#lowage <- targage - hw; highage <- targage + hw
min_age <- abinbeg-(abw/2); max_age <- abinend #+(abw/2)
nbins <- length(abinage)
yfit <- matrix(NA, nrow=nbins, ncol=maxreps)

# arrays for sample number tracking
ndec <- matrix(0, ncol=nbins, nrow=ns)
ndec_tot <- rep(0, nbins)
#xspan <- rep(0, ntarg)
ninwin <- matrix(0, ncol=nbins, nrow=ns)
ninwin_tot <- rep(0, nbins)



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
        influx[1:nsamp,ii] <- indata$zt # 
        nsamples[ii] <- nsamp
    }
  }
}
nsites <- ii

# number of sites with data
nsites

# trim samples to age range
influx[age >= abinend+abw/2] <- NA
age[age >= abinend+abw/2] <- NA

# censor abs(influx) values > 10
influx[abs(influx) >= 10] <- NA
age[abs(influx) >= 10] <- NA


# count number of sites that contributed to each fitted value
ptm <- proc.time()
for (i in 1:nbins) {
  
  for (j in 1:nsites) {
    for (k in 1:nsamples[j]) {
      if (!is.na(age[k,j])) {
        ii <- as.integer(ceiling((age[k,j]-abinbeg-(abw/2.0))/abw))+1
        #print (c(i,j,k,ii))
        if (ii > 0 && ii <= nbins) {ndec[j,ii] = 1}
        if (age[k,j] >= abinage[i]-(abw/2) && age[k,j] <= abinage[i]+(abw/2)) {ninwin[j,i] = 1}
      }
    }
  }
  ndec_tot[i] <- sum(ndec[,i])
  ninwin_tot[i] <- sum(ninwin[,i])
  # xspan[i] <- agemax - agemin
}
proc.time() - ptm
head(cbind(1950 - abinage, pltYearCE, abinage,ndec_tot,ninwin_tot))
tail(cbind(1950 - abinage, pltYearCE, abinage,ndec_tot,ninwin_tot))

ptm <- proc.time()
# 1. reshape matrices into vectors 
x <- as.vector(age)
y <- as.vector(influx)
lfdata <- data.frame(x,y)
lfdata <- na.omit(lfdata)
lfdata <- lfdata[lfdata$x >= min_age & lfdata$x < max_age, ]
x <- lfdata$x; y <- lfdata$y

# average influx for each age bin
binnum <- as.integer(ceiling((x-abinbeg-(abw/2.0))/abw))+1
binave <- tapply(y, binnum, mean)
binaveage <- tapply(x, binnum, mean)
bin_fit_all <- binave

head(cbind(1950 - abinage, pltYearCE, abinage, bin_fit_all))
tail(cbind(1950 - abinage, pltYearCE, abinage, bin_fit_all))

# Bootstrap samples

# Step 1 -- Set up to plot individual replications
if (plotout == "pdf") {pdf(file=paste(curvecsvpath,pdffile,sep=""))}
plot(NULL, ylim=ylim2, xlim=xlim, ylab="Z-scores Tranformed Influx", xlab=xlab, 
     cex.sub=0.8, sub=curvename, type="n")
axis(side = 1, at = seq(xmin, xmax, by = xminortick), labels = FALSE, tcl = -.25) 
axis(side = 1, at = seq(xmin, xmax, by = xminortick*5), labels = FALSE, tcl = -.5) 

# for debugging step plots
# plot the bin averages -- note pltYearCE offset to center the plot steps
# points(1950 - x, y, pch=16, cex=0.5, col=rgb(0.5,0.5,0.5,0.70))
# lines(bin_fit_all ~ pltYearCE, type="s", col="red", lwd=2)
# points(1950 - abinage, bin_fit_all, col="blue", pch=16, cex=0.5)

set.seed(10) # do this to get the same sequence of random samples for each run

# Step 2 -- Do the bootstrap iterations, and plot each age-bin curve
ptm <- proc.time() # time the loop
for (i in 1:nreps) {
  #print(i)
  randsitenum <- sample(seq(1:nsites), nsites, replace=TRUE)
  # print(head(randsitenum))
  x <- as.vector(age[,randsitenum])
  y <- as.vector(influx[,randsitenum])
  lfdata <- data.frame(x,y)
  lfdata <- na.omit(lfdata)
  lfdata <- lfdata[lfdata$x >= min_age & lfdata$x < max_age, ]
  x <- lfdata$x; y <- lfdata$y
  
  binnum <- as.integer(ceiling((x-abinbeg-(abw/2.0))/abw))+1
  binave <- tapply(y, binnum, mean)
  binaveage <- tapply(x, binnum, mean)
  bin_fit <- binave
  yfit[,i] <- bin_fit
  lines(pltYearCE, yfit[,i], lwd=2, type="s", col=rgb(0.5,0.5,0.5,0.10))
}
proc.time() - ptm # how long?

# Step 3 -- Plot the unresampled (initial) area averages
lines(pltYearCE, bin_fit_all, type="s", lwd=2, col="red")

# Step 4 -- Find and add bootstrap CIs
yfit975 <- apply(yfit, 1, function(x) quantile(x,prob=0.975, na.rm=T))
yfit025 <- apply(yfit, 1, function(x) quantile(x,prob=0.025, na.rm=T))
yfit50 <- apply(yfit, 1, function(x) quantile(x,prob=0.500, na.rm=T))
lines(pltYearCE, yfit975, type="s", lwd=1, col="red")
lines(pltYearCE, yfit025, type="s", lwd=1, col="red")

if (plotout == "pdf") {dev.off()}
curveout <- data.frame(cbind(abinage, 1950-abinage, bin_fit_all, yfit975, yfit025, ndec_tot, ninwin_tot))
colnames(curveout) <- c("age", "YearCE", "bin_ave", "cu95", "cl95", "nsites", "ninwin")
outputfile <- paste(curvecsvpath, curvefile, sep="")
write.table(curveout, outputfile, col.names=TRUE, row.names=FALSE, sep=",")




