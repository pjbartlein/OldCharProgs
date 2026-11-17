# composite curve via block binning
# bootstrap-by-site confidence intervals

# paths for input and output .csv files -- modify as appropriate
datapath <- "e:/Projects/GPWG/GPWGv3/data/v3d/"
outpath <- "e:/Projects/GPWG/GPWGv3/data/v3h/v3h2_blocks_R/"

# regionfile
regionpath <- "e:/Projects/GPWG/GPWGv3/data/v3d/v3h2_regions/"
regionname <- "v3h2_nsa"

# names
queryname <- "v3d"
basename <- "zt2kb"
binname <- "bw20" # prebinning

# number of bootstrap replications
nreps <- 200

# target ages for fitted values
targbeg <- -50
targend <- 2200
targstep <- 50
ylabel <- "bin_50"

# plot output (pdf or screen)
plotout <- "pdf" # "screen" # 

# no changes below here

# presampled/binned files
csvpath <- "e:/Projects/GPWG/GPWGv3/data/v3d/v3d_presamp_csv_R/"
csvname <- paste("_presamp_influx_",basename,"_",binname,".csv", sep="")
csvname

# regionfile
regionfile <- paste(regionpath, regionname, ".csv", sep="")
regionfile

# curve (output) file
curvecsvpath <- paste(datapath,queryname,"_blocks_R/",sep="")
curvecsvpath
# if output folder does not exist, create it
dir.create(file.path(curvecsvpath), showWarnings=FALSE)
curvefile <- paste(regionname,"_locfit_",basename,"_",binname,"_",as.character(hw),"_",
  as.character(nreps), ".csv", sep="")
print(curvefile)

# .pdf plot of bootstrap iterations
if (plotout == "pdf") {
  pdffile <- paste(regionname,"_locfit_",basename,"_",binname,"_",as.character(hw),"_",
    as.character(nreps),".pdf", sep="")
  print(pdffile)
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
targage <- seq(targbeg,targend,targstep)
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
influx[age >= targend+targstep/2] <- NA
age[age >= targend+targstep/2] <- NA

# censor abs(influx) values > 10
influx[abs(influx) >= 10] <- NA
age[abs(influx) >= 10] <- NA


# reshape matrices into vectors 
x <- as.vector(age)
y <- as.vector(influx)
lfdata <- data.frame(x,y)
lfdata <- na.omit(lfdata)
x <- lfdata$x; y <- lfdata$y

# write.table(lfdata, "debug.csv", col.names=TRUE, row.names=FALSE, sep=",")

# initial fit, unresampled (i.e. all) data
#loc01 <- locfit(y ~ lp(x, deg=1, h=hw), maxk=800, family="qrgauss")
#summary(loc01)

binnum <- as.integer(ceiling((x-targbeg-(targstep/2))/targstep))+1
binave <- tapply(y, binnum, mean)
binaveage <- tapply(x, binnum, mean)
binave.fit <- binave

# plot the bin averages
plot(targage-(targstep/2), binave, xlab="age", ylab=ylabel, xlim=c(0,2200), ylim=c(-1.0,1.0), 
  col="red", lwd=2, type="s")


# Bootstrap-by-site confidence intervals

# Step 1 -- Set up to plot individual replications

if (plotout == "pdf") {pdf(file=paste(curvecsvpath,pdffile,sep=""))}
plot(targage-(targstep/2), binave, xlab="age", ylab=ylabel, 
  xlim=c(0,2200), ylim=c(-1.0,1.0), type="n")

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
  
  binnum <- as.integer(ceiling((x-targbeg-(targstep/2))/targstep))+1
  binave <- tapply(y, binnum, mean)
  yfit[,i] <- binave
  
  # note plotting lines is slowww
  lines(targage-(targstep/2), yfit[,i], lwd=2, type="s", col=rgb(0.5,0.5,0.5,0.10))
  if (i %% 10 == 0) {print(i)}
}

# Step 3 -- Plot the unresampled (initial) fit
lines(targage-(targstep/2), binave.fit, type="s", lwd=2, col="red")

# Step 4 -- Find and add bootstrap CIs
yfit95 <- apply(yfit, 1, function(x) quantile(x,prob=0.975, na.rm=T))
yfit05 <- apply(yfit, 1, function(x) quantile(x,prob=0.025, na.rm=T))
lines(targage-(targstep/2), yfit95, type="s", lwd=1, col="red")
lines(targage-(targstep/2), yfit05, type="s", lwd=1, col="red")
if (plotout == "pdf") {dev.off()}

# output
curveout <- data.frame(cbind(targage, binave.fit, yfit95, yfit05))
colnames(curveout) <- c("age", "locfit", "cu95", "cl95")
outputfile <- paste(curvecsvpath, curvefile, sep="")
write.table(curveout, outputfile, col.names=TRUE, row.names=FALSE, sep=",")
warnings()