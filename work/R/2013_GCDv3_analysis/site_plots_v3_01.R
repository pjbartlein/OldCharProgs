# simple plots of each site

#csvpath <- "e:/Projects/GPWG/GPWGv3/data/v3h/v3h_trans_csv_fortran/"
#presamppath <- "/Projects/GPWG/GPWGv3/data/v3h/v3h_presamp_csv/"
#pdfpath <- "e:/Projects/GPWG/GPWGv3/data/v3h/v3h_sites_pdf/"
#sitenamefile <- "/Projects/GPWG/GPWGv3/data/v3h/v3h_regions/v3h_sites.csv"
#queryname <- "v3h"

csvfile <- "/Projects/GPWG/GPWGv3/data/v3d-Jul2013/v3d_trans/"
presamppath <- "/Projects/GPWG/GPWGv3/data/v3d-Jul2013/v3d_presamp/"
pdfpath <- "e:/Projects/GPWG/GPWGv3/data/v3d-Jul2013/v3d_sites_pdf/"
sitenamefile <- "/Projects/GPWG/GPWGv3/data/v3d-Jul2013/v3d_regions/v3d_all.csv"
queryname <- "v3d"

sitenames <- read.csv(sitenamefile, stringsAsFactors=FALSE)
sitenames$SiteName <- as.character(sitenames$Sitename)
nsites <- max(sitenames[,1])


for (j in seq(1,nsites)) {	
	sitenum <- j
	print(sitenum)
	siteidchar <- as.character(sitenum)
	print(siteidchar)
	if (sitenum >= 1) siteid <- paste("000", siteidchar, sep="")
	if (sitenum >= 10) siteid <- paste("00", siteidchar, sep="")
	if (sitenum >= 100) siteid <- paste("0", siteidchar, sep="")
  if (sitenum >= 1000) siteid <- paste("", siteidchar, sep="")
	print(siteid)
	csvfile <- paste(csvpath, siteid, "_trans_influx_zt21k.csv", sep="")
	print(csvfile)
  presampfile <- paste(presamppath, siteid, "_presamp_influx_zt21k_sm10.csv", sep="")
  print (presampfile)
	pdfname <- paste(pdfpath, siteid, "_site.pdf", sep="")
	print(pdfname)

	if (file.exists(csvfile)) {
		sitedata <- read.csv(csvfile, header=TRUE, row.names=NULL) 
    presampdata <- read.csv(presampfile, header=TRUE, row.names=NULL) 
		print(csvfile) 
		jj <- which(sitenames[,1] == j)
		print(sitenames$Sitename[jj])

		pdf(file=pdfname ,width=8.5, height=11.0)
		#windows(width=8.5, height=11.0)
		
		# plot raw data
	  par(fig=c(0,1,0.8,1.0), mar=c(2.1, 4.1, 2.1, 4.1))
		plot(sitedata[,6] ~ sitedata[,4], xlim=c(20000,-100), axes=F, mgp=c(2,0,0),
			main=paste(queryname,": ",siteidchar, " ", sitenames$SiteName[jj]), font.main=1, lab=c(10,5,5), 
			ylab="quant", xlab="Age", cex.lab=0.8, pch=16, cex=0.5, type="o")
		axis(1, cex.axis=0.8, xaxp=c(0,22000,22)); axis(2, cex.axis=0.8); axis(4, cex.axis=0.8)
    
	  # plot transformed data
	  par(new=T, fig=c(0,1,0.6,0.8), mar=c(2.1, 4.1, 2.1, 4.1))
	  plot(sitedata[,8] ~ sitedata[,4], xlim=c(20000,-100), axes=F, mgp=c(2,0,0), lab=c(10,5,5),
	    ylab="tall", xlab="Age", cex.lab=0.8, pch=16, cex=0.5, type="o")
	  axis(1, cex.axis=0.8, xaxp=c(0,22000,22)); axis(2, cex.axis=0.8); axis(4, cex.axis=0.8)

	  # z-scores 		
	  par(new=T, fig=c(0,1,0.4,0.6), mar=c(2.1, 4.1, 2.1, 4.1))
	  plot(sitedata[,9] ~ sitedata[,4], xlim=c(20000,-100), axes=F, mgp=c(2,0,0), lab=c(10,5,5), 
	    ylab="Z-Score", xlab="Age", cex.lab=0.8, pch=16, cex=0.5, type="o")
	  axis(1, cex.axis=0.8, xaxp=c(0,22000,22)); axis(2, cex.axis=0.8); axis(4, cex.axis=0.8)
	  
	  if (length(presampdata[,1]) > 0) {
  	  # presample z-scores   	
  	  par(new=T, fig=c(0,1,0.2,0.4), mar=c(2.1, 4.1, 2.1, 4.1))
  	  plot(presampdata[,3] ~ presampdata[,2], xlim=c(20000,-100), axes=F, mgp=c(2,0,0), lab=c(10,5,5), 
  	    ylab="Presample Z-Score", xlab="Age", cex.lab=0.8, pch=16, cex=0.5, type="o")
  	  axis(1, cex.axis=0.8, xaxp=c(0,22000,22)); axis(2, cex.axis=0.8); axis(4, cex.axis=0.8)
	}	  
	dev.off()
	}	
}



