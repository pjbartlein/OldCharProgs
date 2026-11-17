workpath <- "e:\\Projects\\GPWG\\work\\global\\figs\\paper03\\SIFig_model_and_ice\\"
djffile <- "tmp2m_djf_zave_anm_2000_lwr.csv"
jjafile <- "tmp2m_jja_zave_anm_2000_lwr.csv"

outpath <- workpath
outname <- "djf_jja_corr.csv"

djf <- read.csv(paste(workpath, djffile, sep=""))
jja <- read.csv(paste(workpath, jjafile, sep=""))

nvars <- length(jja)
cross_corr <- matrix(nrow=nvars, ncol=1)
corr_pvalue <- matrix(nrow=nvars, ncol=1)

for (i in 1:nvars) {
    cross_corr[i] <- cor(djf[,i], jja[,i], use="na.or.complete")
    corr_pvalue[i] <- cor.test(djf[,i], jja[,i], use="na.or.complete")$p.value
}

rownames(cross_corr) <- names(jja[,1:nvars])

write.table(cbind(cross_corr,corr_pvalue), file=paste(outpath, outname, sep=""), row.names=TRUE, col.names=FALSE, sep=" ")



prefile <- "pre_land_zave_anm_2000_lwr.csv"
tmp2mfile <- "tmp2m_land_zave_anm_2000_lwr.csv"

outpath <- workpath
outname <- "tmp2m_pre_corr.csv"

pre <- read.csv(paste(workpath, prefile, sep=""))
tmp2m <- read.csv(paste(workpath, tmp2mfile, sep=""))

nvars <- length(pre)
cross_corr <- matrix(nrow=nvars, ncol=1)
corr_pvalue <- matrix(nrow=nvars, ncol=1)

for (i in 1:nvars) {
    cross_corr[i] <- cor(pre[,i], tmp2m[,i], use="na.or.complete")
    corr_pvalue[i] <- cor.test(pre[,i], tmp2m[,i], use="na.or.complete")$p.value
}

rownames(cross_corr) <- names(pre[,1:nvars])

write.table(cbind(cross_corr,corr_pvalue), file=paste(outpath, outname, sep=""), row.names=TRUE, col.names=FALSE, sep=" ")