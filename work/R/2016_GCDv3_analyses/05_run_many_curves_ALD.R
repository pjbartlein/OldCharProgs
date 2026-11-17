# Run curves
rm(list=ls())

source("C:/Users/daniau/Dropbox/GCDv3plus/R_code/04_smooth_curve05_jm.R")

runcurves(region="NH", latmin=0, latmax=90, lonmin= -180, lonmax= 180, sitelist="v3plus_all_noTagamaucia.csv", rootpath="C:/Users/daniau/Dropbox/GCDv3plus/")
runcurves("WNA")
runcurves("SEUR")
runcurves("Patag")
runcurves("SEAust")


if(region=="ENA"){sites2 <- sites[sites$Lat>25 & sites$Lat<52 & sites$Lon>(-100) & sites$Lon<(-50),]}
if(region=="BOR"){sites2 <- sites[sites$Lat>52 & sites$Lat<68 & sites$Lon>(-100) & sites$Lon<(-50),]}
if(region=="NEUR"){sites2 <- sites[sites$Lat>52 & sites$Lat<70 & sites$Lon>(-12) & sites$Lon<31,]}
if(region=="SEUR"){sites2 <- sites[sites$Lat>36.7 & sites$Lat<52 & sites$Lon>(-12) & sites$Lon<29,]}
if(region=="Patag"){sites2 <- sites[sites$Lat>(-56) & sites$Lat<(-32) & sites$Lon>(-80) & sites$Lon<(-60),]}
if(region=="SEAust"){sites2 <- sites[sites$Lat>(-44) & sites$Lat<(-21) & sites$Lon>133 & sites$Lon<155,]}
if(region=="WNA"){sites2 <- sites[sites$Lat>30 & sites$Lat<50 & sites$Lon>(-125) & sites$Lon<(-100),]}
if(region=="Globe"){sites2 <- sites}
