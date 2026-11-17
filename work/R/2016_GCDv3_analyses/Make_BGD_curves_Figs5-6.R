# BGD Marlon et al. 2016
######################################## 

setwd("~/data/projects/gcdv3/R/")
#pfSiteSel()
data(paleofiredata)
data(paleofiresites)

############# FIGURE 5 #################

a1 <- pfSiteSel(lat>0)
all=pfTransform(a1,method=c("MinMax","Box-Cox","MinMax","Z-Score"),BasePeriod=c(200,1000),QuantType="INFL")
COMP2=pfCompositeLF(all, tarAge=seq(-60,1000,5), hw=10, nboot=199)
plot(COMP2, main="NH Biomass Burning [base: -10,200]", xlim=c(1000,-100))
result1 <- data.frame(COMP2$Result)
result1$Year <- 1950-(result1$AGE)
names(result1) <- c("YearBP","CharcoalIndex","BootMean","CI05","CI95","Year")
write.csv(result1, paste0("composites/BBNH_pre5_10hw_1k.csv"), row.names=F)

a2 <- pfSiteSel(lat<0)
all=pfTransform(a2,method=c("MinMax","Box-Cox","MinMax","Z-Score"),BasePeriod=c(200,1000),QuantType="INFL")
COMP2=pfCompositeLF(all, tarAge=seq(-60,1000,5), hw=10, nboot=199)
plot(COMP2, main="SH Biomass Burning [base: -10,200]", xlim=c(1000,-100))
result1 <- data.frame(COMP2$Result)
result1$Year <- 1950-(result1$AGE)
names(result1) <- c("YearBP","CharcoalIndex","BootMean","CI05","CI95","Year")
write.csv(result1, paste0("composites/BBSH_pre5_10hw_1k.csv"), row.names=F)

a3 <- pfSiteSel()
all=pfTransform(a3,method=c("MinMax","Box-Cox","MinMax","Z-Score"),BasePeriod=c(200,1000),QuantType="INFL")
COMP2=pfCompositeLF(all, tarAge=seq(-60,1000,5), hw=10, nboot=199)
plot(COMP2, main="GL Biomass Burning [base: -10,200]", xlim=c(1000,-100))
result1 <- data.frame(COMP2$Result)
result1$Year <- 1950-(result1$AGE)
names(result1) <- c("YearBP","CharcoalIndex","BootMean","CI05","CI95","Year")
write.csv(result1, paste0("composites/BBGL_pre5_10hw_1k.csv"), row.names=F)

#wider smoothing window
a1a <- pfSiteSel(lat>0)
all=pfTransform(a1a,method=c("MinMax","Box-Cox","MinMax","Z-Score"),BasePeriod=c(200,1000),QuantType="INFL")
COMP2=pfCompositeLF(all, tarAge=seq(-60,1000,10), hw=20, nboot=199)
plot(COMP2, main="NH Biomass Burning [base: 200-1k]", xlim=c(1000,-100))
result1 <- data.frame(COMP2$Result)
result1$Year <- 1950-(result1$AGE)
names(result1) <- c("YearBP","CharcoalIndex","BootMean","CI05","CI95","Year")
write.csv(result1, paste0("composites/BBNH_pre10_20hw_1k.csv"), row.names=F)

a2a <- pfSiteSel(lat<0)
all=pfTransform(a2a,method=c("MinMax","Box-Cox","MinMax","Z-Score"),BasePeriod=c(200,1000),QuantType="INFL")
COMP2=pfCompositeLF(all, tarAge=seq(-60,1000,10), hw=20, nboot=199)
plot(COMP2, main="SH Biomass Burning [base: 200-1k]", xlim=c(1000,-100))
result1 <- data.frame(COMP2$Result)
result1$Year <- 1950-(result1$AGE)
names(result1) <- c("YearBP","CharcoalIndex","BootMean","CI05","CI95","Year")
write.csv(result1, paste0("composites/BBSH_pre10_20hw_1k.csv"), row.names=F)

a3a <- pfSiteSel()
all=pfTransform(a3a,method=c("MinMax","Box-Cox","MinMax","Z-Score"),BasePeriod=c(200,1000),QuantType="INFL")
COMP2=pfCompositeLF(all, tarAge=seq(-60,1000,10), hw=20, nboot=199)
plot(COMP2, main="GL Biomass Burning [base: 200-1k]", xlim=c(1000,-100))
result1 <- data.frame(COMP2$Result)
result1$Year <- 1950-(result1$AGE)
names(result1) <- c("YearBP","CharcoalIndex","BootMean","CI05","CI95","Year")
write.csv(result1, paste0("composites/BBGL_pre10_20hw_1k.csv"), row.names=F)




############# June 27, 2015 FIGURE 7 #################
write.csv(paleofiresites, file="paleofiresitesV3_R.csv")
# used GCDv3_MapData_Fig1.csv to identify v2 sites, then made a new file with just GCDv2 sites

# ---- global v2
v2 <- read.csv("paleofiresitesV2.csv")
#m <- pfSiteSel() #n=736, v3
m <- pfSiteSel(site_name %in% as.character(v2$site_name)) #n=665
title="GCDv2 Global"
all=pfTransform(m,method=c("MinMax","Box-Cox","MinMax","Z-Score"),BasePeriod=c(200,22000),QuantType="INFL")
COMP2=pfCompositeLF(all, tarAge=seq(-60,22000,20), hw=500, nboot=199)
plot(COMP2, main=paste0(title, " Biomass Burning [base: 200,22000]"), xlim=c(22000,-100), ylim=c(-5,5))
save(all, file="GBB200-22000b_500hw_data_v2.Rda")
save(COMP2, file="GBB200-22000b_500hw_curve_v2.Rda")
result1 <- data.frame(COMP2$Result)
result1$Year <- 1950-(result1$AGE)
names(result1) <- c("YearBP","CharcoalIndex","BootMean","CI05","CI95","Year")
write.csv(result1, paste0("composites/GBB_pre20_500hw_22k_v2.csv"), row.names=F)

m <- pfSiteSel() #n=736, v3
#m <- pfSiteSel(site_name %in% as.character(v2$site_name)) #n=665
title="GCDv3 Global"
all=pfTransform(m,method=c("MinMax","Box-Cox","MinMax","Z-Score"),BasePeriod=c(200,22000),QuantType="INFL")
COMP2=pfCompositeLF(all, tarAge=seq(-60,22000,20), hw=500, nboot=199)
plot(COMP2, main=paste0(title, " Biomass Burning [base: 200,22000]"), xlim=c(22000,-100), ylim=c(-5,5))
save(all, file="GBB200-22000b_500hw_data_v3.Rda")
save(COMP2, file="GBB200-22000b_500hw_curve_v3.Rda")
result1 <- data.frame(COMP2$Result)
result1$Year <- 1950-(result1$AGE)
names(result1) <- c("YearBP","CharcoalIndex","BootMean","CI05","CI95","Year")
write.csv(result1, paste0("composites/GBB_pre20_500hw_22k_v3.csv"), row.names=F)



# ---- southern hemisphere (extratropics)
#v2 <- read.csv("paleofiresitesV2.csv")
m <- pfSiteSel(lat< -29.99) #n=186, v3
#m <- pfSiteSel(lat< -29.99, site_name %in% as.character(v2$site_name)) #n=140
title="GCDv3 Southern Extratropics"
all=pfTransform(m,method=c("MinMax","Box-Cox","MinMax","Z-Score"),BasePeriod=c(200,22000),QuantType="INFL")
COMP2=pfCompositeLF(all, tarAge=seq(-60,22000,20), hw=500, nboot=199)
plot(COMP2, main=paste0(title, " Biomass Burning [base: 200,22000]"), xlim=c(22000,-100), ylim=c(-5,5))
save(all, file="SHBB200-22000b_500hw_data_v3.Rda")
save(COMP2, file="SHBB200-22000b_500hw_curve_v3.Rda")
result1 <- data.frame(COMP2$Result)
result1$Year <- 1950-(result1$AGE)
names(result1) <- c("YearBP","CharcoalIndex","BootMean","CI05","CI95","Year")
write.csv(result1, paste0("composites/SHBB_pre20_500hw_22k_v3.csv"), row.names=F)

# ---- tropics
#v2 <- read.csv("paleofiresitesV2.csv")
#m <- pfSiteSel(lat> -29.99 & lat < 30.01) #n=200, v3
m <- pfSiteSel(lat> -29.99 & lat < 30.01, site_name %in% as.character(v2$site_name)) #n=665
title="GCDv2 Tropics"
all=pfTransform(m,method=c("MinMax","Box-Cox","MinMax","Z-Score"),BasePeriod=c(200,22000),QuantType="INFL")
COMP2=pfCompositeLF(all, tarAge=seq(-60,22000,20), hw=500, nboot=199)
plot(COMP2, main=paste0(title, " Biomass Burning [base: 200,22000]"), xlim=c(22000,-100), ylim=c(-5,5))
save(all, file="TRBB200-22000b_500hw_data_v2.Rda")
save(COMP2, file="TRBB200-22000b_500hw_curve_v2.Rda")
result1 <- data.frame(COMP2$Result)
result1$Year <- 1950-(result1$AGE)
names(result1) <- c("YearBP","CharcoalIndex","BootMean","CI05","CI95","Year")
write.csv(result1, paste0("composites/TRBB_pre20_500hw_22k_v2.csv"), row.names=F)

# ---- Northern hemisphere (extratropics)
v2 <- read.csv("paleofiresitesV2.csv")
#m <- pfSiteSel(lat> 29.99) #n=736, v3
m <- pfSiteSel(lat> 29.99, site_name %in% as.character(v2$site_name)) #n=665
title="GCDv2 Northern Extratropics"
all=pfTransform(m,method=c("MinMax","Box-Cox","MinMax","Z-Score"),BasePeriod=c(200,22000),QuantType="INFL")
COMP2=pfCompositeLF(all, tarAge=seq(-60,22000,20), hw=500, nboot=199)
plot(COMP2, main=paste0(title, " Biomass Burning [base: 200,22000]"), xlim=c(22000,-100), ylim=c(-5,5))
save(all, file="NHBB200-22000b_500hw_data_v2.Rda")
save(COMP2, file="NHBB200-22000b_500hw_curve_v2.Rda")
result1 <- data.frame(COMP2$Result)
result1$Year <- 1950-(result1$AGE)
names(result1) <- c("YearBP","CharcoalIndex","BootMean","CI05","CI95","Year")
write.csv(result1, paste0("composites/NHBB_pre20_500hw_22k_v2.csv"), row.names=F)

m <- pfSiteSel(lat> 29.99) #n=736, v3
#m <- pfSiteSel(lat> 29.99, site_name %in% as.character(v2$site_name)) #n=665
title="GCDv3 Northern Extratropics"
all=pfTransform(m,method=c("MinMax","Box-Cox","MinMax","Z-Score"),BasePeriod=c(200,22000),QuantType="INFL")
COMP2=pfCompositeLF(all, tarAge=seq(-60,22000,20), hw=500, nboot=199)
plot(COMP2, main=paste0(title, " Biomass Burning [base: 200,22000]"), xlim=c(22000,-100), ylim=c(-5,5))
save(all, file="NHBB200-22000b_500hw_data_v3.Rda")
save(COMP2, file="NHBB200-22000b_500hw_curve_v3.Rda")
result1 <- data.frame(COMP2$Result)
result1$Year <- 1950-(result1$AGE)
names(result1) <- c("YearBP","CharcoalIndex","BootMean","CI05","CI95","Year")
write.csv(result1, paste0("composites/NHBB_pre20_500hw_22k_v3.csv"), row.names=F)




