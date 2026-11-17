# pca missing data
library(pcaMethods)

pcapath="e:\\Projects\\GPWG\\analysis\\output\\v2n\\v2n_pcafiles\\"
aname="pca02"
window="0500"
pnmiss="50"

subsetfile <- paste(pcapath, aname, "_", window, "_", pnmiss, "_subset.csv", sep="")
locfile <- paste(pcapath, aname, "_", window, "_", pnmiss, ".dat", sep="")

zin <- read.csv(subsetfile)
#zin <- read.csv(file.choose()) # NA's = blanks
mapdat <-  read.table(locfile)
#mapdat <- read.table(file.choose())

lat <- mapdat[,1]
lon <- mapdat[,2]
age <- zin[,1]
zmat <- data.matrix(zin)[,2:(length(zin)-2)]

zd <- prep(zmat, scale="none", center=FALSE)

# PCA
pcamethod="ppca"
ncomp <- 10

resPCA <- pca(zd, method=pcamethod, center=FALSE, nPcs=ncomp)
resPCA
#plotPcs(resPCA,pcs=1:nP(resPCA), type="loadings")

loadpca <- loadings(resPCA)
scores <- scores(resPCA)
loadnames <- rownames(loadpca)

meanhighloadvars <- matrix(nrow=resPCA@nObs, ncol=resPCA@nPcs)
loadcor <- matrix(nrow=1, ncol=resPCA@nPcs)
signcor <- matrix(nrow=1, ncol=resPCA@nPcs)

for (k in 1:resPCA@nPcs) {
    #k=1
    nhighloadvars <- 10
    highloadvars <- matrix(nrow=resPCA@nObs, ncol=nhighload)
    rankload <- rank(loadpca[,k])
    rankload <- length(rankload)-rankload+1 #reverse order
    for (j in 1:nhighloadvars) {
        #j=1
        highloadvars[,j] <- zmat[,loadnames[rankload==j]]
    }
    meanhighloadvars[,k] <- rowMeans(highloadvars, na.rm=TRUE)
    
    plot(age,meanhighloadvars[,k])
    plot(age,scores[,k]/resPCA@sDev[k])
    plot(meanhighloadvars[,k],scores[,k]/resPCA@sDev[k])
    
    loadcor[1,k] <- cor(meanhighloadvars[,k],scores[,k], use="pairwise.complete.obs")
}
loadcor
signcor <- sign(loadcor)
signcor

scoresout <- cbind(age,scores)
loadout <- cbind(lat,lon,loadpca)

ncompchar <- paste ("0", as.character(ncomp), sep="")
if (ncomp >= 10) ncompchar <- as.character(ncomp)
scoresfile <- paste(pcapath, aname, "_", window, "_", pnmiss, "_", ncompchar,"_scores.csv", sep="")
loadfile <- paste(pcapath, aname, "_", window, "_", pnmiss, "_", ncompchar,"_load.dat", sep="")
write.table(scoresout, file=scoresfile, row.names=FALSE, col.names=TRUE, sep=",")
write.table(loadout, file=loadfile, row.names=FALSE, col.names=FALSE, sep=" ")

#write.table(scoresout, file=file.choose(), row.names=FALSE, col.names=TRUE, sep=",")
#write.table(loadout, file=file.choose(),row.names=FALSE, col.names=FALSE, sep=" ")

