# BGD Marlon et al. 2016
######################################## 

library(GCD)
library(paleofire)
data(paleofiredata)
data(paleofiresites)

# attempt to reproduce Fig. 5 "Globe" curve

# get data and transform
a3 <- pfSiteSel() # Globe
all=pfTransform(a3,method=c("MinMax","Box-Cox","MinMax","Z-Score"),BasePeriod=c(200,1000),QuantType="INFL")
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
  write.csv(transdata, paste(trans_csvpath, siteid, "_trans_influx_pf.csv"), row.names=F)
}

COMP2=pfCompositeLF(all, tarAge=seq(-60,1000,5), hw=10, nboot=199, pseudodata = FALSE)
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
  write.csv(prebinneddata, paste(prebin_csvpath, siteid, "_prebin_influx_pf.csv"), row.names=F)
}

# why does this say "base:  -10, 200)?
plot(COMP2, main="GL Biomass Burning [base: -10,200]", xlim=c(1000,-100))
result1 <- data.frame(COMP2$Result)
result1$Year <- 1950-(result1$AGE)
curves_csvpath <- "e:/Projects/GPWG/GPWGv3/GCDv3Data/v3i_Rpackage/v3i_curves/"
names(result1) <- c("YearBP","CharcoalIndex","BootMean","CI05","CI95","Year")
write.csv(result1, paste(curves_csvpath, "BBGL_pre5_10hw_1k.csv", sep=""), row.names=F)


COMP2b=pfCompositeLF(all, tarAge=seq(-60,1000,5), hw=20, nboot=199, pseudodata=FALSE)
attributes(COMP2b)

# why does this say "base:  -10, 200)?
plot(COMP2b, main="GL Biomass Burning [base: -10,200]", xlim=c(1000,-100))
result2 <- data.frame(COMP2b$Result)
result2$Year <- 1950-(result2$AGE)
curves_csvpath <- "e:/Projects/GPWG/GPWGv3/GCDv3Data/v3i_Rpackage/v3i_curves/"
names(result2) <- c("YearBP","CharcoalIndex","BootMean","CI05","CI95","Year")
write.csv(result2, paste(curves_csvpath, "BBGL_pre5_20hw_1k.csv", sep=""), row.names=F)

# read Fortran version composite curves
Fortran_csvpath="e:/Projects/GPWG/GPWGv3/GCDv3Data/v3i_Fortran/v3i_curves/"
fort_hw10 <- read.csv(paste(Fortran_csvpath, "paleofiresitesV3_R_boot_influx_zt1k_sm02_0010.csv", sep=""))
fort_hw10$Year <- 1950-(fort_hw10$targ_age)

plot(COMP2, main="GL Biomass Burning [base: -10,200] -- hw10", xlim=c(1000,-100))
lines(fort_hw10$targ_age, fort_hw10$paleofiresitesV3_R_influx_zt1k_sm02_0010, col="blue", lwd=2)
lines(fort_hw10$targ_age, fort_hw10$yfit_boot_cu95, col="blue")
lines(fort_hw10$targ_age, fort_hw10$yfit_boot_cl95, col="blue")

Fortran_csvpath="e:/Projects/GPWG/GPWGv3/GCDv3Data/v3i_Fortran/v3i_curves/"
fort_hw20 <- read.csv(paste(Fortran_csvpath, "paleofiresitesV3_R_boot_influx_zt1k_sm02_0020.csv", sep=""))
fort_hw20$Year <- 1950-(fort_hw20$targ_age)

plot(COMP2b, main="GL Biomass Burning [base: -10,200 -- hw20]", xlim=c(1000,-100))
lines(fort_hw20$targ_age, fort_hw20$paleofiresitesV3_R_influx_zt1k_sm02_0020, col="blue", lwd=2)
lines(fort_hw20$targ_age, fort_hw20$yfit_boot_cu95, col="blue")
lines(fort_hw20$targ_age, fort_hw20$yfit_boot_cl95, col="blue")

diff_10hw <- fort_hw10$paleofiresitesV3_R_influx_zt1k_sm02_0010 - result1$CharcoalIndex 
diff_20hw <- fort_hw20$paleofiresitesV3_R_influx_zt1k_sm02_0020 - result2$CharcoalIndex 
plot(result1$YearBP, diff_10hw, xlim=c(1050, -100), xaxp  = c(1000, -100, 11), type="l", col="magenta")
lines(result1$YearBP, diff_20hw, col="purple")

cbind(result1$Year, result1$YearBP, diff_10hw, diff_20hw)

# R Scripts
RScripts_csvpath="e:/Projects/GPWG/GPWGv3/GCDv3Data/v3i_RScripts/v3i_curves/"
Rscripts_hw20 <- read.csv(paste(RScripts_csvpath, "paleofiresitesV3_R.csv_locfit_zt1k_bw5_20_199.csv", sep=""))
Rscripts_hw20$Year <- 1950-(Rscripts_hw20$age)

# reproduce Fig. 5
plot(result1$Year, result1$CharcoalIndex, xlim=c(1000,2000), xaxp  = c(1000, 2000, 10), ylim=c(-1,1), col="black", type="l", lwd=2, 
     main="Fig. 5")
lines(result2$Year, result2$CharcoalIndex, col="red", lwd=2)
lines(result2$Year, result2$CI05, col="red", lwd=1)
lines(result2$Year, result2$CI95, col="red", lwd=1)
lines(fort_hw20$Year, fort_hw20$paleofiresitesV3_R_influx_zt1k_sm02_0020, col="blue", lwd=2)
lines(fort_hw20$Year, fort_hw20$yfit_boot_cu95, col="blue")
lines(fort_hw20$Year, fort_hw20$yfit_boot_cl95, col="blue")
lines(Rscripts_hw20$Year, Rscripts_hw20$locfit, col="purple", lwd=2)
lines(Rscripts_hw20$Year, Rscripts_hw20$cu95, col="purple")
lines(Rscripts_hw20$Year, Rscripts_hw20$cl95, col="purple")
