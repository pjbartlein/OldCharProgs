#nclr <- 12
#plotclr <- rgb(
#	c(0.600,0.800,1.000,1.000,1.000,1.000,0.887,0.607,0.376,0.194,0.000,0.000),
#	c(0.050,0.067,0.250,0.501,0.740,0.924,1.000,0.918,0.792,0.630,0.167,0.133),
#	c(0.000,0.000,0.000,0.144,0.376,0.694,1.000,1.000,1.000,1.000,1.000,0.800),
#	maxColorValue=1)
	
nclr <- 5
plotclr <- rgb(
	c(1.000,1.000,0.800,0.376,0.000),
	c(0.000,0.501,0.800,0.792,0.000),
	c(0.000,0.144,0.800,1.000,1.000),
	maxColorValue=1)

csvpath <- "\\Projects\\GPWG\\data\\gcd_v2\\csv_files\\trans_influx2_zts_csv_v2h\\"
pdfpath <- "\\Projects\\GPWG\\data\\plots\\site_plots\\"
sitenamefile <- "\\Projects\\GPWG\\data\\gcd_v2\\files\\misc\\v2h_sites_elev.csv"
meansfile <- "\\Projects\\GPWG\\data\\gcd_v2\\files\\slices\\v2h_means.csv"
anomsfile <- "\\Projects\\GPWG\\data\\gcd_v2\\files\\slices\\v2h_anoms.csv"


sitenames <- read.csv(sitenamefile, stringsAsFactors=FALSE)
sitenames$SiteName <- as.character(sitenames$Sitename)
nsites <- max(sitelist$ID_SI)

meansin <- read.csv(meansfile)
means <- as.matrix(meansin[,4:47])
anomsin <- read.csv(anomsfile)
anoms <- as.matrix(anomsin[,4:47])

j <- 1
for (j in seq(1,nsites)) {
#for (j in seq(1,1)) {	
	sitenum <- j
	print(sitenum)
	siteidchar <- as.character(sitenum)
	print(siteidchar)
	if (sitenum >= 1) siteid <- paste("000", siteidchar, sep="")
	if (sitenum >= 10) siteid <- paste("00", siteidchar, sep="")
	if (sitenum >= 100) siteid <- paste("0", siteidchar, sep="")
	print(siteid)
	csvfile <- paste(csvpath, siteid, "_trans_influx2_zts.csv", sep="")
	print(csvfile)
	pdfname <- paste(pdfpath, siteid, "_site.pdf", sep="")
	print(pdfname)

	if (file.exists(csvfile)) {
		sitedata <- read.csv(csvfile) 
		print(csvfile) 
		jj <- which(sitenames$ID_SI == j)
		print(sitenames$Sitename[jj])

		pdf(file=pdfname ,width=8.5, height=11.0)
		#windows(width=8.5, height=11.0)
		
		# plot raw data
		par(fig=c(0,1,.75,1), mar=c(2.1, 4.1, 2.1, 4.1))
		plot(sitedata$influx2 ~ sitedata$est_age, xlim=c(20000,-100), axes=F, mgp=c(2,0,0),
			main=paste(siteidchar, " -- ", sitenames$SiteName[jj]), font.main=1, lab=c(10,5,5), 
			ylab=paste("Influx2", sitenames$PREF_UNITS[jj]), xlab="Age", cex.lab=0.8, pch=16, cex=0.5, type="o")
		axis(1, cex.axis=0.8, xaxp=c(0,22000,22)); axis(2, cex.axis=0.8); axis(4, cex.axis=0.8)
		
		# plot transformed data
		par(new=T, fig=c(0,1,.5,.75), mar=c(2.1, 4.1, 2.1, 4.1))
		plot(sitedata$tall ~ sitedata$est_age, xlim=c(20000,-100), axes=F, mgp=c(2,0,0), lab=c(10,5,5),
			ylab="tall", xlab="Age", cex.lab=0.8, pch=16, cex=0.5, type="o")
		axis(1, cex.axis=0.8, xaxp=c(0,22000,22)); axis(2, cex.axis=0.8); axis(4, cex.axis=0.8)
		
		# ages
		age <- seq(21000,500, by=-500)
		age[43]=450.0
		age[44]=0.0		

		# z-scores and means			
		par(new=T, fig=c(0,1,.25,.5), mar=c(2.1, 4.1, 2.1, 4.1))
		plot(sitedata$ztpi ~ sitedata$est_age, xlim=c(20000,-100), axes=F, mgp=c(2,0,0), lab=c(10,5,5), 
			ylab="Z-Score", xlab="Age", cex.lab=0.8, pch=16, cex=0.5, type="l", col="gray")
		axis(1, cex.axis=0.8, xaxp=c(0,22000,22)); axis(2, cex.axis=0.8); axis(4, cex.axis=0.8)
		
		#colornum <- cut(means[jj,], c(-999.0,-5.0,-4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,999.0), labels=FALSE)
		colornum <- cut(means[jj,], c(-999.0,-1.15,-0.375,0.375,1.15,999.0), labels=FALSE)
		colcode <- plotclr[nclr - colornum + 1]
		points(sitedata$ztpi ~ sitedata$est_age, pch=16, cex=0.6, col="gray")
		if (max(means[jj,], na.rm=T) > -9999.0) {
		for (k in c(seq(1,42, by=1),44) ) {
			points(x=age[k],y=means[jj,k], col=colcode[k], pch=16, cex=0.9)
			points(x=age[k],y=means[jj,k], cex=0.9) }}
		
		#  anomalies
		#mean1k <- NA
		#tempage <- sitedata$est_age[sitedata$est_age < 700.0 & sitedata$est_age > 200.0 & sitedata$ztpi != -9999.0]
		#tempdat <- sitedata$ztpi[sitedata$est_age < 700.0 & sitedata$est_age > 200.0 & sitedata$ztpi != -9999.0]
		#if (length(tempdat) != 0) mean1k <- mean(tempdat)		
		#print(mean1k)
		anm <- rep(0.0, length(sitedata$tall))
		#if (!is.na(mean1k)) anm <- sitedata$ztpi - mean1k
		if (means[jj,43] != -9999.0) anm <- sitedata$ztpi - means[jj,43]

		par(new=T, fig=c(0,1,0,.25), mar=c(2.1, 4.1, 2.1, 4.1))
		plot(anm ~ sitedata$est_age, xlim=c(20000,-100), axes=F, mgp=c(2,0,0), lab=c(10,5,5), 
			ylab="Anomaly", xlab="Age", cex.lab=0.8, pch=16, cex=0.5, type="l", col="gray")
		#plot(anm ~ sitedata$est_age, xlim=c(20000,-100), ylim=c(-3,3), axes=F, mgp=c(2,0,0), lab=c(10,5,5), type="n",
		#	ylab="Anomaly", xlab="Age")
		points(anm ~ sitedata$est_age, pch=16, cex=0.6, col="gray")


		axis(1, cex.axis=0.8, xaxp=c(0,22000,22)); axis(2, cex.axis=0.8); axis(4, cex.axis=0.8)
		#colornum <- cut(sitedata$anm, c(-999.0,-5.0,-4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,999.0), labels=FALSE)
		#colornum <- cut(means[jj,], c(-999.0,-1.15,-0.375,0.375,1.15,999.0), labels=FALSE)
		#colcode <- plotclr[nclr - colornum + 1]
					
		#colornum <- cut(anoms[jj,], c(-999.0,-5.0,-4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,999.0), labels=FALSE)
		colornum <- cut(anoms[jj,], c(-999.0,-1.15,-0.375,0.375,1.15,999.0), labels=FALSE)
		colcode <- plotclr[nclr - colornum + 1]
		
		plotchar <- 24
		if (max(anoms[jj,], na.rm=T) > -9999.0) {
		for (k in c(seq(1,42, by=1),44) ) {
			if (anoms[jj,k] > -9999.0) {
			yval <- anoms[jj,k]
			xval <- age[k]
			if (yval <= -0.375) points(x=xval, y=yval, bg=colcode[k], pch=25, cex=0.8)
			if (yval >= 0.375) points(x=xval, y=yval, bg=colcode[k], pch=24, cex=0.8)
			if (yval > -0.375 & yval < 0.375) points(x=xval, y=yval, bg=colcode[k], pch=22, cex=0.8) 
			} } }	

	dev.off()
	}
}



