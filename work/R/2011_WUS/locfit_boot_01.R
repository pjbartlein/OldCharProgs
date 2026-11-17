# locfit
library(locfit)

workpath <-"e:\\Projects\\IMPD\\2010\\data\\"
firescarfile <- "scar_w104.dat"
infile <- paste(workpath, firescarfile, sep="")
outname <- "locfit_w104_01.csv"
outfile <- paste(workpath, outname, sep="")

# locfit params
nn_param <- .2
#acri_param <- "GCV"

# matrix for fitted values
nfityrs <- 800
maxreps <- 1000
yfit <- matrix(nrow=nfityrs,ncol=maxreps)
newYear <- data.frame(x=seq(1201, 2000,  1))

# read the data
fsread <- read.table(infile)
nsites <- length(fsread)-1
sitenum <- seq(1,nsites,1) # integer sitenums, not database site numbers

# detach Year column
Yearin <- fsread[,1]
fsin <- fsread[,1:nsites+1]
nyrs <- length(Yearin)

# data for last iteration, omit observations with pscar NaN's
NReport <- apply(fsin,1,FUN=function(x) length(which(x >= 0)) )
NScar <- apply(fsin,1,FUN=function(x) length(which(x > 0)) )
pscar <- NScar/NReport
Year <- Yearin[pscar[]!="NaN"]
PropScar <- pscar[pscar[]!="NaN"]
plot(PropScar ~ Year, col="gray", pch=16, cex=0.5)

nreps <- 1000
set.seed(10) # dp this to get the same sequence of samples
plot(PropScar ~ Year, col="black", pch=16, cex=0.5)
for (i in 1:nreps) {
    print(i)

    # random sample of site numbers
    #fs <- fsin  # no sampling, just the data as it comes in
    randsitenum <- sample(sitenum,nsites,replace=T)
    fs <- fsin[,randsitenum]  # new dataset
    
    # omit observations with pscar Nan's
    nreport <- apply(fs,1,FUN=function(x) length(which(x >= 0)) )
    nscar <- apply(fs,1,FUN=function(x) length(which(x > 0)) )
    pscar <- nscar/nreport
    x <- Yearin[pscar[]!="NaN"]
    y <- pscar[pscar[]!="NaN"]
    
    # locfit
    loc01 <- locfit(y ~ lp(x, deg=2, nn=nn_param), family="binomial")
    #print(summary(loc01))
    pred01 <- predict(loc01, newdata=newYear, se.fit=F)
    lines(newYear$x, pred01, col=rgb(.5,.5,.5,0.10), lwd=3)
    yfit[,i] <- pred01
}

# last fit -- original (unresampled) data
x <- Year
y <- PropScar
loc01 <- locfit(y ~ lp(x, deg=2, nn=nn_param), family="binomial")
#print(summary(loc01))
pred01 <- predict(loc01, newdata=newYear, se.fit=F)
#plot(PropScar ~ Year, col="gray", pch=16, cex=0.5)
lines(newYear$x, pred01, col="red", lwd=3)

yfit95 <- apply(yfit,1,function(x) quantile(x,prob=0.95,na.rm=T))
yfit05 <- apply(yfit,1,function(x) quantile(x,prob=0.05,na.rm=T))
lines(newYear$x, yfit95, col="red", lwd=1)
lines(newYear$x, yfit05, col="red", lwd=1)

# output
targetyr <- newYear$x
out.locfit <- cbind(targetyr,pred01,yfit05,yfit95)
write.table(out.locfit, file=outfile, row.names=FALSE, col.names=TRUE, sep=",")
