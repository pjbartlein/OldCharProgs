pcapath="e:\\Projects\\GPWG\\analysis\\output\\v2n\\v2n_pcafiles\\"
aname="z22-0nrs"
window="0500"
pnmiss="50"
pcamethod <- "ppca"
ncomp="04"

scorespltfile <- paste(pcapath, aname, "_", window, "_", pnmiss, "_", pcamethod, "_", ncomp, "_scores.csv", sep="")
plottitle <- paste(aname, "_", window, "_", pnmiss, "_", pcamethod, "_", ncomp, "_scores.csv", sep="")
scoresplt <- read.csv(scorespltfile)
plottitle

plot(cbind(scoresplt["PC1"],scoresplt[c("globe","nhemis","shemis","nextr","sextr","ntrop",
    "strop","hemis_diff","extr_diff","trop_diff","nlatdiff","slatdiff")]),main=plottitle)            

plot(cbind(scoresplt["PC2"],scoresplt[c("globe","nhemis","shemis","nextr","sextr","ntrop",
    "strop","hemis_diff","extr_diff","trop_diff","nlatdiff","slatdiff")]),main=plottitle)

plot(cbind(scoresplt["PC3"],scoresplt[c("globe","nhemis","shemis","nextr","sextr","ntrop",
    "strop","hemis_diff","extr_diff","trop_diff","nlatdiff","slatdiff")]),main=plottitle)
            

plot(cbind(scoresplt["PC4"],scoresplt[c("globe","nhemis","shemis","nextr","sextr","ntrop",
    "strop","hemis_diff","extr_diff","trop_diff","nlatdiff","slatdiff")]),main=plottitle)

cmat <- matrix(nrow=12, ncol=4)
cmat[ 1,] <- cbind(cor(scoresplt["globe"],scoresplt["PC1"]),     cor(scoresplt["globe"],scoresplt["PC2"]),     cor(scoresplt["globe"],scoresplt["PC3"]),     cor(scoresplt["globe"],scoresplt["PC4"]))
cmat[ 2,] <- cbind(cor(scoresplt["nhemis"],scoresplt["PC1"]),    cor(scoresplt["nhemis"],scoresplt["PC2"]),    cor(scoresplt["nhemis"],scoresplt["PC3"]),    cor(scoresplt["nhemis"],scoresplt["PC4"]))
cmat[ 3,] <- cbind(cor(scoresplt["shemis"],scoresplt["PC1"]),    cor(scoresplt["shemis"],scoresplt["PC2"]),    cor(scoresplt["shemis"],scoresplt["PC3"]),    cor(scoresplt["shemis"],scoresplt["PC4"]))
cmat[ 4,] <- cbind(cor(scoresplt["nextr"],scoresplt["PC1"]),     cor(scoresplt["nextr"],scoresplt["PC2"]),     cor(scoresplt["nextr"],scoresplt["PC3"]),     cor(scoresplt["nextr"],scoresplt["PC4"]))
cmat[ 5,] <- cbind(cor(scoresplt["sextr"],scoresplt["PC1"]),     cor(scoresplt["sextr"],scoresplt["PC2"]),     cor(scoresplt["sextr"],scoresplt["PC3"]),     cor(scoresplt["sextr"],scoresplt["PC4"]))
cmat[ 6,] <- cbind(cor(scoresplt["ntrop"],scoresplt["PC1"]),     cor(scoresplt["ntrop"],scoresplt["PC2"]),     cor(scoresplt["ntrop"],scoresplt["PC3"]),     cor(scoresplt["ntrop"],scoresplt["PC4"]))
cmat[ 7,] <- cbind(cor(scoresplt["strop"],scoresplt["PC1"]),     cor(scoresplt["strop"],scoresplt["PC2"]),     cor(scoresplt["strop"],scoresplt["PC3"]),     cor(scoresplt["strop"],scoresplt["PC4"]))
cmat[ 8,] <- cbind(cor(scoresplt["hemis_diff"],scoresplt["PC1"]),cor(scoresplt["hemis_diff"],scoresplt["PC2"]),cor(scoresplt["hemis_diff"],scoresplt["PC3"]),cor(scoresplt["hemis_diff"],scoresplt["PC4"]))
cmat[ 9,] <- cbind(cor(scoresplt["extr_diff"],scoresplt["PC1"]), cor(scoresplt["extr_diff"],scoresplt["PC2"]), cor(scoresplt["extr_diff"],scoresplt["PC3"]), cor(scoresplt["extr_diff"],scoresplt["PC4"]))
cmat[10,] <- cbind(cor(scoresplt["trop_diff"],scoresplt["PC1"]), cor(scoresplt["trop_diff"],scoresplt["PC2"]), cor(scoresplt["trop_diff"],scoresplt["PC3"]), cor(scoresplt["trop_diff"],scoresplt["PC4"]))
cmat[11,] <- cbind(cor(scoresplt["nlatdiff"],scoresplt["PC1"]),  cor(scoresplt["nlatdiff"],scoresplt["PC2"]),  cor(scoresplt["nlatdiff"],scoresplt["PC3"]),  cor(scoresplt["nlatdiff"],scoresplt["PC4"]))
cmat[12,] <- cbind(cor(scoresplt["slatdiff"],scoresplt["PC1"]),  cor(scoresplt["slatdiff"],scoresplt["PC2"]),  cor(scoresplt["slatdiff"],scoresplt["PC3"]),  cor(scoresplt["slatdiff"],scoresplt["PC4"]))
rownames(cmat) <- c("globe","nhemis","shemis","nextr","sextr","ntrop",
    "strop","hemis_diff","extr_diff","trop_diff","nlatdiff","slatdiff")
colnames(cmat) <- c("PC1", "PC2", "PC3", "PC4")
cmat

