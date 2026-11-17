# composite curve via block binning
# bootstrap-by-site confidence intervals

# paths for input and output .csv files -- modify as appropriate
datapath <- "e:/Projects/GPWG/work/feedback/data/ice_core/"
outpath <- datapath

# number of bootstrap replications
nreps <- 10

# target ages for fitted values
targbeg <- -50
targend <- 2200
targstep <- 50
ylabel <- "bin_50"

dataname <- "law2006_ch4"

# plot output (pdf or screen)
plotout <- "pdf" # "screen" # 

# no changes below here

inputname <- paste(dataname, ".csv", sep="")
inputname

# curve (output) file
curvefile <- paste(dataname,"_",as.character(targstep),"_", as.character(nreps), ".csv", sep="")
print(curvefile)

# .pdf plot of bootstrap iterations
if (plotout == "pdf") {
  pdffile <- paste(dataname,"_",as.character(targstep),"_", as.character(nreps),".pdf", sep="")
  print(pdffile)
}

# matrices for data and fitted values
targage <- seq(targbeg,targend,targstep)
length(targage)
maxreps <- 1000
yfit <- matrix(NA, nrow=length(targage), ncol=maxreps)

# read CH4 data
inputfile <- paste(datapath, inputname, sep="")
print(inputfile)
indata <- read.csv(inputfile)
head(indata)

age <- indata$Age_yrBP1950 # 
CH4 <- indata$CH4_ppb_NOAA04_Scale # 

# trim samples to age range
CH4[age >= targend+targstep/2] <- NA
age[age >= targend+targstep/2] <- NA

# reshape matrices into vectors 
x <- as.vector(age)
y <- as.vector(CH4)
lfdata <- data.frame(x,y)
lfdata <- na.omit(lfdata)
x <- lfdata$x; y <- lfdata$y
one <- rep(1,length(x))

binnum <- as.integer(ceiling((x-targbeg-(targstep/2))/targstep))+1
binave <- tapply(y, binnum, mean)
binaveage <- tapply(x, binnum, mean)
bincount <- tapply(one, binnum, sum)
binave.fit <- binave

# bin numbers of each bin that had an average (or a single) value
binsub <- as.numeric(unlist(dimnames(binave)))


ageout <- matrix(NA, ncol=1, nrow=length(targage))
binaveout <- matrix(NA, ncol=1, nrow=length(targage))
bincountout <- matrix(NA, ncol=1, nrow=length(targage))

ageout <- targage
binaveout[binsub] <- binave[binsub]
bincountout[binsub] <- bincount[binsub]
binout <- data.frame(cbind(ageout,binaveout,bincountout))
names(binout) <- c("targ_age", "CH4", "count")

# plot the bin averages
pdf(file=paste(outpath, "d13C_binned.pdf", sep=""))
plot(binout$targ_age-(targstep/2), binout$CH4, xlab="age", ylab=ylabel, xlim=c(0,2200), 
  col="black", lwd=2, type="s")
dev.off()

# output\
outputfile <- paste(outpath, curvefile, sep="")
write.table(binout, outputfile, col.names=TRUE, row.names=FALSE, sep=",")
warnings()