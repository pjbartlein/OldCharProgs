# correlations and pseudo-loadings

pcapath <- "e:\\Projects\\GPWG\\analysis\\output\\v2n\\v2n_pcafiles\\"
aname <- "z22-0"
hw <- "0500"
pnmiss <- "00"

curvespath <- "e:\\Projects\\GPWG\\analysis\\output\\v2n\\v2n_pcafiles\\"
curvesname <- "z22-0nrs_0500_50_ppca_04_scores.csv"

outpath <- "e:\\Projects\\GPWG\\analysis\\output\\v2n\\v2n_pcafiles\\"
outname <- "pc_scores"

subsetfile <- paste(pcapath, aname, "_", hw, "_", pnmiss, "_subset.csv", sep="")
curvesfile <- paste(curvespath, curvesname,  sep="")
locfile <- paste(pcapath, aname, "_", hw, "_", pnmiss, ".dat", sep="")
outfile <- paste(outpath, aname, "nrs_", hw, "_", pnmiss, "_", outname, "_corr.dat", sep="")

zin <- read.csv(subsetfile)
curvesin <- read.csv(curvesfile)
#zin <- read.csv(file.choose()) # NA's = blanks
mapdat <-  read.table(locfile)
#mapdat <- read.table(file.choose())

lat <- mapdat[,1]
lon <- mapdat[,2]
age <- zin[,1]

ncurves <- 4
zmat <- data.matrix(zin)[,2:(length(zin)-2)]
curvesmat <- data.matrix(curvesin)[,2:(ncurves+1)]

nvars <- length(zin)-3
cross_corr <- matrix(nrow=nvars, ncol=ncurves)

for (i in 1:nvars) {
    for (j in 1:ncurves) {
#    i <- 1
#    j <- 1
        cross_corr[i,j] <- cor(zmat[,i], curvesmat[,j], use="na.or.complete")
    }
}

colnames(cross_corr) <- names(curvesin[,2:(ncurves+1)])
corr_out <- cbind(lat,lon,cross_corr)

write.table(corr_out, file=outfile, row.names=FALSE, col.names=FALSE, sep=" ")

