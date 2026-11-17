# composite curve via the locfit package
# bootstrap-by-site confidence intervals

library(locfit)

# paths for input and output .csv files -- modify as appropriate
datapath <- "e:/Projects/GPWG/GPWGv3/data/v3d_Jul2013/"
outpath <- "e:/Projects/GPWG/GPWGv3/data/v3h/v3h_curves_R_csv/"

# regionfile
regionpath <- "e:/Projects/GPWG/GPWGv3/data/v3h/v3h_regions/"
regionname <- "v3h_nsa" 

# names
queryname <- "v3h"
basename <- "zt21k"
binname <- "bw20" # 

# locfit window-width parameter
hw <- 200 # bandwidth (smoothing parameter)
# number of bootstrap replications
nreps <- 10

# target ages for fitted values
targbeg <- -60
targend <- 22000
targstep <- 20

# plot output (pdf or screen)
plotout <- "pdf" # "screen" # 

# no changes below here

# presampled/binned files
#csvpath <- paste(datapath,queryname,"_presample_csv/", sep="")
csvpath <- "e:/Projects/GPWG/GPWGv3/data/v3d_Jul2013/v3d_presample_csv_R/"
csvname <- paste("_presamp_influx_",basename,"_",binname,".csv", sep="")

# regionfile
regionfile <- paste(regionpath, regionname, ".csv", sep="")

# curve (output) file
curvecsvpath <- paste(datapath,queryname,"_curves_R_csv/",sep="")
# if output folder does not exist, create it
dir.create(file.path(datapath, paste(queryname,"_curves_R_csv/",sep="")), showWarnings=FALSE)
curvefile <- paste(regionname,"_locfit_",basename,"_",binname,"_",as.character(hw),"_",
  as.character(nreps),".csv", sep="")

# .pdf plot of bootstrap iterations
if (plotout == "pdf") {
  pdffile <- paste(regionname,"_locfit_",basename,"_",binname,"_",as.character(hw),"_",
    as.character(nreps),".pdf", sep="")
}

# read the list of sites
sites <- read.csv(regionfile)
head(sites)

ns <- length(sites[,1]) #length(sites$ID_SITE)
ns

# matrices for data and fitted values
maxrecs <- 2000
age <- matrix(NA, ncol=ns, nrow=maxrecs)
influx <- matrix(NA, ncol=ns, nrow=maxrecs)
targage <- data.frame(x=seq(targbeg,targend,targstep))
length(targage$x)
maxreps <- 1000
yfit <- matrix(NA, nrow=length(targage$x), ncol=maxreps)

# read and store the presample (binned) files as matrices of ages and influx values
ii <- 0
for (i in 1:ns) {
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
    }
  }
}
nsites <- ii
nsites
#head(cbind(age,influx))


# trim samples to age range
influx[age >= targend+hw] <- NA
age[age >= targend+hw] <- NA

# censor abs(influx) values > 10
influx[abs(influx) >= 10] <- NA
age[abs(influx) >= 10] <- NA

# locfit

# reshape matrices into vectors 
x <- as.vector(age)
y <- as.vector(influx)
lfdata <- data.frame(x,y)
lfdata <- na.omit(lfdata)
x <- lfdata$x; y <- lfdata$y

# write.table(lfdata, "debug.csv", col.names=TRUE, row.names=FALSE, sep=",")

# initial fit, unresampled (i.e. all) data
loc01 <- locfit(y ~ lp(x, deg=1, h=hw), maxk=800, family="qrgauss")
summary(loc01)

# get and plot fitted values
pred01 <- predict(loc01, newdata=targage, se.fit=TRUE)
locfit <- data.frame(targage$x, pred01$fit)
fitname <- paste("locfit_",as.character(hw), sep="")
colnames(locfit) <- c("age", fitname)
head(locfit)
plot(locfit[,1], locfit[,2], xlab="age", ylab=paste(regionname,fitname,sep=" "),
  type="l", lwd=2, col="red", ylim=c(-2,2))


# Bootstrap-by-site confidence intervals

# Step 1 -- Set up to plot individual replications
if (plotout == "pdf") {pdf(file=paste(curvecsvpath,pdffile,sep=""))}
plot(x, y, xlab="age", ylab=paste(regionname,fitname,sep=" "), 
  ylim=c(-0.8,0.8), xlim=c(0,22000), type="n")

# Step 2 -- Do the bootstrap iterations, and plot each composite curve
set.seed(10) # do this to get the same sequence of random samples for each run

for (i in 1:nreps) {
  #print(i)
  randsitenum <- sample(seq(1:nsites), nsites, replace=TRUE)
  # print(head(randsitenum))
  x <- as.vector(age[,randsitenum])
  y <- as.vector(influx[,randsitenum])
  lfdata <- data.frame(x,y)
  lfdata <- na.omit(lfdata)
  x <- lfdata$x; y <- lfdata$y
  locboot <- locfit(y ~ lp(x, deg=1, h=hw), maxk=800, maxit=20, family="qrgauss")
  predboot <- predict(locboot, newdata=targage, se.fit=TRUE)
  yfit[,i] <- predboot$fit
  # note plotting lines is slowww
  lines(targage$x, yfit[,i], lwd=2, col=rgb(0.5,0.5,0.5,0.10))
  if (i %% 10 == 0) {print(i)}
}

# Step 3 -- Plot the unresampled (initial) fit
locfit <- data.frame(targage$x, pred01$fit)
fitname <- paste("locfit_",as.character(hw), sep="")
colnames(locfit) <- c("age", fitname)
lines(locfit[,1], locfit[,2], lwd=2, col="red")

# Step 4 -- Find and add bootstrap CIs
yfit95 <- apply(yfit, 1, function(x) quantile(x,prob=0.975, na.rm=T))
yfit05 <- apply(yfit, 1, function(x) quantile(x,prob=0.025, na.rm=T))
lines(targage$x, yfit95, lwd=1, col="red")
lines(targage$x, yfit05, lwd=1, col="red")
if (plotout == "pdf") {dev.off()}

# output
curveout <- data.frame(cbind(targage$x, pred01$fit, yfit95, yfit05))
colnames(curveout) <- c("age", "locfit", "cu95", "cl95")
outputfile <- paste(curvecsvpath, curvefile, sep="")
write.table(curveout, outputfile, col.names=TRUE, row.names=FALSE, sep=",")

