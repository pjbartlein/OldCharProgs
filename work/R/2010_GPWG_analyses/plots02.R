# correlations between PC scores and regional composite curves

pcapath <-"e:\\Projects\\GPWG\\analysis\\output\\v2n\\v2n_pcafiles\\"
curvesname <- "v2n_all_nsa_boot_influx_zt21k_sm10_0500.csv"
aname <- "z22-0nrs"
hw <- "0500"
pnmiss <- "50"
pcamethod <- "ppca"
ncomp <- 4
ncompchar <- paste ("0", as.character(ncomp), sep="")
if (ncomp >= 10) ncompchar <- as.character(ncomp)

scoresfile <- paste(pcapath,aname,"_",hw,"_",pnmiss,"_",pcamethod,"_",ncompchar,"_scores.csv",sep="")
curvesfile <- paste(pcapath,curvesname, sep="")
outfile <- paste(pcapath,aname,"_",hw,"_",pnmiss,"_",pcamethod,"_",ncompchar,"_curvepca.csv",sep="")

scoresin <- read.csv(scoresfile)
curvesin <- read.csv(curvesfile)
sc <- cbind(scoresin,curvesin)

plottitle <- paste(aname, "_", window, "_", pnmiss, "_", pcamethod, "_", ncompchar, "_scores.csv", sep="")
plottitle

plot(cbind(sc["PC1"],sc[c("globe","nhemis","shemis","nextr","sextr","ntrop",
    "strop","hemis_diff","extr_diff","trop_diff","nlatdiff","slatdiff")]),main=plottitle)            

plot(cbind(sc["PC2"],sc[c("globe","nhemis","shemis","nextr","sextr","ntrop",
    "strop","hemis_diff","extr_diff","trop_diff","nlatdiff","slatdiff")]),main=plottitle)

plot(cbind(sc["PC3"],sc[c("globe","nhemis","shemis","nextr","sextr","ntrop",
    "strop","hemis_diff","extr_diff","trop_diff","nlatdiff","slatdiff")]),main=plottitle)
            

plot(cbind(sc["PC4"],sc[c("globe","nhemis","shemis","nextr","sextr","ntrop",
    "strop","hemis_diff","extr_diff","trop_diff","nlatdiff","slatdiff")]),main=plottitle)

cmat <- matrix(nrow=12, ncol=4)
cmat[ 1,] <- cbind(cor(sc["globe"],sc["PC1"]),     cor(sc["globe"],sc["PC2"]),     cor(sc["globe"],sc["PC3"]),     cor(sc["globe"],sc["PC4"]))
cmat[ 2,] <- cbind(cor(sc["nhemis"],sc["PC1"]),    cor(sc["nhemis"],sc["PC2"]),    cor(sc["nhemis"],sc["PC3"]),    cor(sc["nhemis"],sc["PC4"]))
cmat[ 3,] <- cbind(cor(sc["shemis"],sc["PC1"]),    cor(sc["shemis"],sc["PC2"]),    cor(sc["shemis"],sc["PC3"]),    cor(sc["shemis"],sc["PC4"]))
cmat[ 4,] <- cbind(cor(sc["nextr"],sc["PC1"]),     cor(sc["nextr"],sc["PC2"]),     cor(sc["nextr"],sc["PC3"]),     cor(sc["nextr"],sc["PC4"]))
cmat[ 5,] <- cbind(cor(sc["sextr"],sc["PC1"]),     cor(sc["sextr"],sc["PC2"]),     cor(sc["sextr"],sc["PC3"]),     cor(sc["sextr"],sc["PC4"]))
cmat[ 6,] <- cbind(cor(sc["ntrop"],sc["PC1"]),     cor(sc["ntrop"],sc["PC2"]),     cor(sc["ntrop"],sc["PC3"]),     cor(sc["ntrop"],sc["PC4"]))
cmat[ 7,] <- cbind(cor(sc["strop"],sc["PC1"]),     cor(sc["strop"],sc["PC2"]),     cor(sc["strop"],sc["PC3"]),     cor(sc["strop"],sc["PC4"]))
cmat[ 8,] <- cbind(cor(sc["hemis_diff"],sc["PC1"]),cor(sc["hemis_diff"],sc["PC2"]),cor(sc["hemis_diff"],sc["PC3"]),cor(sc["hemis_diff"],sc["PC4"]))
cmat[ 9,] <- cbind(cor(sc["extr_diff"],sc["PC1"]), cor(sc["extr_diff"],sc["PC2"]), cor(sc["extr_diff"],sc["PC3"]), cor(sc["extr_diff"],sc["PC4"]))
cmat[10,] <- cbind(cor(sc["trop_diff"],sc["PC1"]), cor(sc["trop_diff"],sc["PC2"]), cor(sc["trop_diff"],sc["PC3"]), cor(sc["trop_diff"],sc["PC4"]))
cmat[11,] <- cbind(cor(sc["nlatdiff"],sc["PC1"]),  cor(sc["nlatdiff"],sc["PC2"]),  cor(sc["nlatdiff"],sc["PC3"]),  cor(sc["nlatdiff"],sc["PC4"]))
cmat[12,] <- cbind(cor(sc["slatdiff"],sc["PC1"]),  cor(sc["slatdiff"],sc["PC2"]),  cor(sc["slatdiff"],sc["PC3"]),  cor(sc["slatdiff"],sc["PC4"]))
rownames(cmat) <- c("globe","nhemis","shemis","nextr","sextr","ntrop",
    "strop","hemis_diff","extr_diff","trop_diff","nlatdiff","slatdiff")
colnames(cmat) <- c("PC1", "PC2", "PC3", "PC4")
cmat

