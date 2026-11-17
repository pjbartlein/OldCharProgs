
csvpath <- "E:\\Projects\\igbpfire\\charsum\\data\\current_hs2\\sites\\site_csv\\"
pdfpath <- "E:\\Projects\\igbpfire\\charsum\\data\\current_hs2\\sites\\transplot\\"

j=1
for (j in seq(1,700)) {	
	sitenum <- j
	#print(sitenum)
	siteidchar <- as.character(sitenum)
	#print(siteidchar)

	if (sitenum >= 1) siteid <- paste("000", siteidchar, "_data", sep="")
	if (sitenum >= 10) siteid <- paste("00", siteidchar, "_data", sep="")
	if (sitenum >= 100) siteid <- paste("0", siteidchar, "_data", sep="")
	#print(siteid)
	csvfile <- paste(csvpath, siteid, ".csv", sep="")
	#print(csvfile)
	pdfname <- paste(pdfpath, siteid, ".pdf", sep="")
	#print(pdfname)

	if (file.exists(csvfile)) {
		sitedata <- read.csv(csvfile) 		
		
		bfit <- boxcox.fit(sitedata$quant, lambda=seq(-2,2,by=0.2), lambda2=TRUE)
		lambda1 <- bfit$lambda[1]
		lambda2 <- bfit$lambda[2]
		y <- sitedata$quant+lambda2
		if (lambda1 == 0.0) quant3 <- log(y) else quant3 <- (y^lambda1 - 1)/lambda1


		pdf(file=pdfname ,width=11.0, height=8.5)
#		windows(width=11.0, height=8.5)
		
		par(fig=c(0,1/3,0,1), mar=c(4.1, 4.1, 2.1, 2.1))
		plot(sitedata$est_age ~ sitedata$quant, ylim=c(22000,-100), axes=F, mgp=c(2,0,0), 
			xlab=paste("quant -- site", sitenum), ylab="age", pch=16, cex=0.8, type="o")
		axis(1); axis(2); axis(3)

		par(new=T, fig=c(1/3,2/3,2/3,1), mar=c(5.1, 4.1, 2.1, 2.1))
		hist(sitedata$quant, xlab="quant", mgp=c(2,1,0), main=NULL)

		par(new=T, fig=c(1/3,2/3,1/3,2/3), mar=c(7.1, 4.1, 0.1, 2.1))
		hist(quant3, xlab=paste("l1=",signif(lambda1,3)," l2= ",signif(lambda2,5)), mgp=c(2,1,0), main=NULL)
		
		par(new=T, fig=c(2/3,1,0,1), mar=c(4.1, 4.1, 2.1, 2.1))
		plot(sitedata$est_age ~ quant3, ylim=c(22000,-100), axes=FALSE, mgp=c(2,0,0), 
			xlab=" transformed quant", ylab="age", pch=16, cex=0.8, type="o")
		axis(1); axis(2); axis(3)

	dev.off()

	}
		
}

	

	
	