# pca missing data
library(pcaMethods)

z500_75 <- read.csv(file.choose()) # NA's = blanks
mapdat_500_75 <- read.table(file.choose())

lat <- mapdat_500_75[,1]
lon <- mapdat_500_75[,2]
age <- z500_75[,1]
zmat <- data.matrix(z500_75)[,2:(length(z500_75)-2)]

z200_50 <- read.csv(file.choose()) # NA's = blanks
mapdat_200_50 <- read.table(file.choose())

lat <- mapdat_200_50[,1]
lon <- mapdat_200_50[,2]
age <- z200_50[,1]
zmat <- data.matrix(z200_50)[,2:(length(z200_50)-2)]


zd <- prep(zmat, scale="none", center=FALSE)

resPPCA <- pca(zd, method="ppca", center=FALSE, nPcs=4)
resPPCA
plotPcs(resPPCA,pcs=1:nP(resPPCA), type="loadings")

resBPCA <- pca(zd, method="bpca", center=FALSE, nPcs=4)
resBPCA
plotPcs(resBPCA,pcs=1:nP(resBPCA), type="loadings")

resSVDI <- pca(zd, method="svdImpute", center=FALSE, nPcs=4)
resSVDI
plotPcs(resSVDI,pcs=1:nP(resSVDI), type="loadings")

resNipals <- pca(zd, method="nipals", center=FALSE, nPcs=4)
resNipals
plotPcs(resNipals,pcs=1:nP(resNipals), type="loadings")

resNLPCA <- pca(zd, method="nlpca", center=FALSE, nPcs=1)
resNLPCA
plotPcs(resNLPCA,pcs=1:nP(resNLPCA), type="loadings")


loadpca <- loadings(resPPCA)
scores <- scores(resPPCA)

scoresout <- cbind(age,scores)
loadout <- cbind(lat,lon,loadpca)
write.table(scoresout, file=file.choose(), row.names=FALSE, col.names=TRUE, sep=",")
write.table(loadout, file=file.choose(),row.names=FALSE, col.names=FALSE, sep=" ")

plot(age,scores(resSVDI)[,1])
plot(age,scores(resSVDI)[,2])
plot(age,scores(resSVDI)[,3])