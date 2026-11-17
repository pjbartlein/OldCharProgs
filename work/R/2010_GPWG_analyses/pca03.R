maxpc <- 10
sDevOut <- matrix(nrow=maxpc,ncol=maxpc)
R2Out <- matrix(nrow=maxpc,ncol=maxpc)
#q2Out <- matrix(nrow=maxpc,ncol=maxpc)
kEstiOut <- matrix(nrow=maxpc,ncol=maxpc)
bestNpcsOut <- matrix(nrow=1, ncol=maxpc)

for (k in 1:maxpc) {
    print(k)
    zd <- prep(zmat, scale="none", center=FALSE)
    resPPCA <- pca(zd, method="ppca", center=FALSE, nPcs=k)
    sDevOut[k,1:k] <- sDev(resPPCA)[1:k]
    R2Out[k,1:k] <- resPPCA@R2[1:k] 
    estifast <- kEstimateFast(zd, method="ppca", em="nrmsep", evalPcs=1:k, allVariables=TRUE,
        verbose=TRUE)
    kEstiOut[k,1:k] <- estifast$eError[1:k]
    bestNpcsOut[1,k] <- estifast$bestNPcs

    rm(resPPCA)
    rm(zd)
}
sDevOut
R2Out
kEstiOut
bestNpcsOut


maxpc <- 12
kEstiOut <- matrix(nrow=maxpc,ncol=maxpc)
bestNpcsOut <- matrix(nrow=1, ncol=maxpc)
for (k in 1:maxpc) {
    print(k)
    zd <- prep(zmat, scale="none", center=FALSE)
    estifast <- kEstimateFast(zd, method="ppca", evalPcs=1:k, allVariables=TRUE)
    kEstiOut[k,1:k] <- estifast$eError[1:k]
    bestNpcsOut[1,k] <- estifast$bestNPcs
}
kEstiOut
bestNpcsOut
