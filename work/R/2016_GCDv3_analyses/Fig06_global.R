# BGD Marlon et al. 2016
######################################## 

library(GCD)
library(paleofire)
data(paleofiredata)
data(paleofiresites)

# attempt to reproduce Fig. 6 "Globe" curve

# get data and transform
m <- pfSiteSel() # Globe
all=pfTransform(m,method=c("MinMax","Box-Cox","MinMax","Z-Score"),BasePeriod=c(200,22000),QuantType="INFL")
attributes(all)
nsites <- dim(all$TransData)[2]
trans_csvpath <- "e:/Projects/GPWG/GPWGv3/GCDv3Data/v3i_Rpackage/v3i_trans_csv/"
for (j in 1:nsites)  {
  
  # get siteid string
  siteidchar <- as.character(paleofiresites$id_site[j])
  if (paleofiresites$id_site[j] >= 1) siteid <- paste("000", siteidchar, sep="")
  if (paleofiresites$id_site[j] >= 10) siteid <- paste("00", siteidchar, sep="")
  if (paleofiresites$id_site[j] >= 100) siteid <- paste("0", siteidchar, sep="")
  if (paleofiresites$id_site[j] >= 1000) siteid <- paste(    siteidchar, sep="")
  
  transdata <- na.omit(data.frame(all$Age[,j], all$IntData[,j], all$TransData[,j]))
  names(transdata) <- c("age", "influx", "zt")
  write.csv(transdata, paste(trans_csvpath, siteid, "_trans_influx_zt21k.csv"), row.names=F)
}

COMP2=pfCompositeLF(all, tarAge=seq(-60,22000,20), hw=500, nboot=199, pseudodata = FALSE)
attributes(COMP2)

prebin_csvpath <- "e:/Projects/GPWG/GPWGv3/GCDv3Data/v3i_Rpackage/v3i_prebin_csv/"
for (j in 1:nsites)  {
  
  # get siteid string
  siteidchar <- as.character(paleofiresites$id_site[j])
  if (paleofiresites$id_site[j] >= 1) siteid <- paste("000", siteidchar, sep="")
  if (paleofiresites$id_site[j] >= 10) siteid <- paste("00", siteidchar, sep="")
  if (paleofiresites$id_site[j] >= 100) siteid <- paste("0", siteidchar, sep="")
  if (paleofiresites$id_site[j] >= 1000) siteid <- paste(    siteidchar, sep="")
  
  prebinneddata <- na.omit(data.frame(COMP2$BinCentres, COMP2$BinnedData[,j]))
  names(prebinneddata) <- c("age", "prebinned")
  write.csv(prebinneddata, paste(prebin_csvpath, siteid, "_prebin_influx_zt21k.csv"), row.names=F)
}

# plot
plot(COMP2, main="GL Biomass Burning", xlim=c(22000,-100))
result1 <- data.frame(COMP2$Result)
curves_csvpath <- "e:/Projects/GPWG/GPWGv3/GCDv3Data/v3i_Rpackage/v3i_curves/"
names(result1) <- c("YearBP","CharcoalIndex","BootMean","CI05","CI95")
write.csv(result1, paste(curves_csvpath, "BBGL_pre20_500hw_zt21k.csv", sep=""), row.names=F)

plot(COMP2, main="GL Biomass Burning [base: -10,200] -- hw10", xlim=c(1000,-100))
lines(fort_hw10$targ_age, fort_hw10$paleofiresitesV3_R_influx_zt1k_sm02_0010, col="blue", lwd=2)
lines(fort_hw10$targ_age, fort_hw10$yfit_boot_cu95, col="blue")
lines(fort_hw10$targ_age, fort_hw10$yfit_boot_cl95, col="blue")

Fortran_csvpath="e:/Projects/GPWG/GPWGv3/GCDv3Data/v3i_Fortran/v3i_curves/"
fort_hw500 <- read.csv(paste(Fortran_csvpath, "paleofiresitesV3_R_boot_influx_zt21k_sm10_0500.csv", sep=""))

diff_500hw <- fort_hw500$paleofiresitesV3_R_influx_zt21k_sm10_0500 - result1$CharcoalIndex 
plot(result1$YearBP, diff_500hw, xlim=c(22000, -100), xaxp  = c(22000, 0, 22), type="l", col="magenta")

cbind(result1$Year, result1$YearBP, diff_10hw, diff_20hw)

# R Scripts
RScripts_csvpath="e:/Projects/GPWG/GPWGv3/GCDv3Data/v3i_RScripts/v3i_curves/"
Rscripts_hw500 <- read.csv(paste(RScripts_csvpath, "paleofiresitesV3_R.csv_locfit_zt21k_bw20_500_199.csv", sep=""))

# reproduce Fig. 6
plot(result1$YearBP, result1$CharcoalIndex, xlim=c(22000,-100), xaxp  = c(22000, 0, 22), ylim=c(-1.5,0.5), col="red", type="l", lwd=2, 
     main="Fig. 6")
lines(result1$YearBP, result1$CI05, col="red", lwd=1)
lines(result1$YearBP, result1$CI95, col="red", lwd=1)
lines(fort_hw500$targ_age, fort_hw500$paleofiresitesV3_R_influx_zt21k_sm10_0500, col="blue", lwd=2)
lines(fort_hw500$targ_age, fort_hw500$yfit_boot_cu95, col="blue")
lines(fort_hw500$targ_age, fort_hw500$yfit_boot_cl95, col="blue")
lines(Rscripts_hw500$age, Rscripts_hw500$locfit, col="purple", lwd=2)
lines(Rscripts_hw500$age, Rscripts_hw500$cu95, col="purple")
lines(Rscripts_hw500$age, Rscripts_hw500$cl95, col="purple")
legend("bottomright", legend=c("Rpackage","Fortran", "Rscripts"), col=c("red", "blue", "purple"), lwd=2)
