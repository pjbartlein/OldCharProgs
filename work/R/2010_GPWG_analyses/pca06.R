# pca missing data
#library(pcaMethods)

pcapath <- "e:\\Projects\\GPWG\\analysis\\output\\v2n\\v2n_pcafiles\\"
aname <- "z15-10"
hw <- "0100"
pnmiss <- "50"

subsetfile <- paste(pcapath, aname, "_", hw, "_", pnmiss, "_subset.csv", sep="")
locfile <- paste(pcapath, aname, "_", hw, "_", pnmiss, ".dat", sep="")

zin <- read.csv(subsetfile)
#zin <- read.csv(file.choose()) # NA's = blanks
mapdat <-  read.table(locfile)
#mapdat <- read.table(file.choose())

lat <- mapdat[,1]
lon <- mapdat[,2]
age <- zin[,1]
zmat <- data.matrix(zin)[,2:(length(zin)-2)]

# PCA
pcamethod="ppca"
ncomp <- 4

zd <- prep(zmat, scale="none", center=FALSE)
set.seed(1)
resPCA <- pca(zd, method=pcamethod, scale="none", center=FALSE, nPcs=ncomp)
resPCA
#plotPcs(resPCA,pcs=1:nP(resPCA), type="loadings")

loadpca <- loadings(resPCA)
scores <- scores(resPCA) # /rep(resPCA@sDev, each=nrow(scores(resPCA)))
scoresout <- cbind(age,scores)
loadout <- cbind(lat,lon,loadpca)
statsout <- cbind(1:resPCA@nPcs,resPCA@sDev,resPCA@R2,resPCA@R2cum)
colnames(statsout) <- c("PC", "sDev", "R2", "R2cum")

ncompchar <- paste ("0", as.character(ncomp), sep="")
if (ncomp >= 10) ncompchar <- as.character(ncomp)
scoresfile <- paste(pcapath,aname,"_",hw,"_",pnmiss,"_",pcamethod,"_",ncompchar,"_scores_ns.csv",sep="")
loadfile <- paste(pcapath,aname,"_",hw,"_",pnmiss,"_",pcamethod,"_",ncompchar,"_load.dat",sep="")
statsfile <- paste(pcapath,aname,"_",hw,"_",pnmiss,"_",pcamethod,"_",ncompchar,"_stats.csv",sep="")

write.table(scoresout, file=scoresfile, row.names=FALSE, col.names=TRUE, sep=",")
write.table(loadout, file=loadfile, row.names=FALSE, col.names=FALSE, sep=" ")
write.table(statsout, file=statsfile, row.names=FALSE, col.names=TRUE, sep=",")
statsout
rm(zd)
rm(resPCA)
rnorm(5)
#rm(zin)
#rm(mapdat)