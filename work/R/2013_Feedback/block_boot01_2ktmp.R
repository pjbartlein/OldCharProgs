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

dataname <- "aaves_tmp_land" # PAGES2k_tave"
ylabel <- paste (dataname,"bin_50", sep=" ")

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

# read data
inputfile <- paste(datapath, inputname, sep="")
print(inputfile)
indata <- read.csv(inputfile)
head(indata)

age <- indata$Age_yrBP1950 # 
tmp2k <- indata$GLOBE # indata$tave # 

# trim samples to age range
tmp2k[age >= targend+targstep/2] <- NA
age[age >= targend+targstep/2] <- NA

# reshape matrices into vectors 
x <- as.vector(age)
y <- as.vector(tmp2k)
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
names(binout) <- c("targage", "Mann_tave", "count")

# plot the bin averages
pdf(file=paste(outpath, "Mann_tave.pdf", sep=""))
plot(binout$targage-(targstep/2), binout$Mann_tave, xlab="age", ylab=ylabel, xlim=c(0,2200), 
  col="green", lwd=2, type="s")
dev.off()

# output\
outputfile <- paste(outpath, curvefile, sep="")
write.table(binout, outputfile, col.names=TRUE, row.names=FALSE, sep=",")
warnings()