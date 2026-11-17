# pca missing data
#library(pcaMethods)

workpath <- "e:\\Projects\\Transient\\work01\\"

tmpname <- "tmp2m_land_zave_anm.csv"
prename <- "pre_land_zave_anm.csv"
pmename <- "p-e_land_zave_anm.csv"

tmp <- read.csv(paste(workpath, tmpname, sep=""))
pre <- read.csv(paste(workpath, prename, sep=""))
pme <- read.csv(paste(workpath, pmename, sep=""))

zin <- cbind(tmp[,50:55], pre[,50:55], pme[,50:55])
age <- tmp[,1]
zmat <- data.matrix(zin)

# PCA via pcaMethods
pcamethod="svd"
ncomp <- 18
zd <- prep(zmat, scale="uv", center=TRUE)
set.seed(1)
resPCA <- pca(zd, method=pcamethod, scale="uv", center=FALSE, nPcs=ncomp)
resPCA
#plotPcs(resPCA,pcs=1:4, type="loadings")
loadpca <- loadings(resPCA)
scores <- scores(resPCA)/rep(resPCA@sDev, each=nrow(scores(resPCA)))

scoresout <- cbind(age,scores)
loadout <- cbind(names(zin),loadpca)
colnames(loadout) <- c("var",colnames(loadpca))
statsout <- cbind(1:resPCA@nPcs,resPCA@sDev,resPCA@R2,resPCA@R2cum)
colnames(statsout) <- c("PC", "sDev", "R2", "R2cum")

# PCA via prcomp
pcamethod="prcomp"
resPCA <- prcomp(zmat, center=TRUE, scale=TRUE )
summary(resPCA)
scores <- predict(resPCA, zmat)/rep(resPCA$sdev, each=nrow(predict(resPCA,zmat)))
loadpca <- resPCA$rotation
unq <- rowSums(loadpca[,1:7]^2)

scoresout <- cbind(age,scores)
loadout <- data.frame(names(zin),loadpca)
colnames(loadout) <- c("var",colnames(loadpca))
#rownames(loadout) <- names(loadpca)
statsout <- cbind(1:length(resPCA$sdev),resPCA$sdev,resPCA$sdev/sum(resPCA$sdev),
    cumsum(resPCA$sdev/sum(resPCA$sdev)))
colnames(statsout) <- c("PC", "sDev", "R2", "R2cum")

# PCA via factanal
pcamethod="fa"
ncomp <- 4
resPCA <- factanal(scale(zmat, center=TRUE, scale=TRUE), control = list(nstart=30,trace=T,lower=0.0001),
    factors=ncomp, scores="regr", rotation="none")
resPCA
scores <- resPCA$scores 
loadpca <- as.matrix(resPCA$loadings[,1:ncomp])

scoresout <- cbind(age,scores)
loadout <- data.frame(names(zin),loadpca)
colnames(loadout) <- c("var",colnames(loadpca))
#rownames(loadout) <- names(loadpca)
statsout <- cbind(1:ncomp,colSums(loadpca^2),colSums(loadpca^2)/length(zin),cumsum(colSums(loadpca^2)/length(zin)))
colnames(statsout) <- c("PC", "sDev", "R2", "R2cum")


# Output
ncompchar <- paste ("0", as.character(ncomp), sep="")
if (ncomp >= 10) ncompchar <- as.character(ncomp)
scoresfile <- paste(workpath,"latzones_",pcamethod,"_",ncompchar,"_scores.csv",sep="")
loadfile <- paste(workpath,"latzones_",pcamethod,"_",ncompchar,"_load.csv",sep="")
statsfile <- paste(workpath,"latzones_",pcamethod,"_",ncompchar,"_stats.csv",sep="")


write.table(scoresout, file=scoresfile, row.names=FALSE, col.names=TRUE, sep=",")
write.table(loadout, file=loadfile, row.names=FALSE, col.names=TRUE, sep=",")
write.table(statsout, file=statsfile, row.names=FALSE, col.names=TRUE, sep=",")
