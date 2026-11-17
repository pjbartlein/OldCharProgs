
csvpath <- "E:\\Projects\\igbpfire\\charsum\\data\\current_hs1\\sites\\sitecsv\\"
pdfpath <- "E:\\Projects\\igbpfire\\charsum\\data\\current_hs1\\sites\\transplot\\"

for (j in seq(1,700)) {	
	sitenum <- j
	#print(sitenum)
	siteidchar <- as.character(sitenum)
	#print(siteidchar)
	if (sitenum >= 1) siteid <- paste("000", siteidchar, sep="")
	if (sitenum >= 10) siteid <- paste("00", siteidchar, sep="")
	if (sitenum >= 100) siteid <- paste("0", siteidchar, sep="")
	#print(siteid)
	csvfile <- paste(csvpath, siteid, ".csv", sep="")
	#print(csvfile)
	pdfname <- paste(pdfpath, siteid, ".pdf", sep="")
	#print(pdfname)

	if (file.exists(csvfile)) {
		sitedata <- read.csv(csvfile) 
		#print(csvfile) 

		char2 <- (sitedata$char-min(sitedata$char))/(max(sitedata$char)-min(sitedata$char))+0.01
		npts <- 101
		y <- char2
		n <- length(y)
		logy <- log(y)
		ydot <- exp(mean(logy))
		lasave <- matrix(1:npts)
		liksave <- matrix(1:npts)
		for (i in 1:npts) {
			la <- -2.0+(i-1)*(4/(npts-1))
			if (la != 0.0) yt <- (y^la-1)/la else yt <- logy*(1+(la*logy)/2*(1+(la*logy)/3*(1+(la*logy)/4)))
			zt <- yt/ydot^(la-1)
			loglik <- -n/2*log(sum((zt - mean(zt))^2 ))
			lasave[i] <- la
			liksave[i] <- loglik
			}
		#print(cbind(lasave,liksave))
		cbind(sitenum,liksave[which.max(liksave)],lasave[which.max(liksave)])
		laopt=as.character(lasave[which.max(liksave)])
		lafit <- lasave[which.max(liksave)]
		if (laopt == 0.0) char3 <- log(char2) else char3 <- (char2^lafit - 1)/lafit


		pdf(file=pdfname ,width=11.0, height=8.5)
		#windows(width=11.0, height=8.5)
		
		par(fig=c(0,1/3,0,1), mar=c(4.1, 4.1, 2.1, 2.1))
		plot(sitedata$age ~ sitedata$char, ylim=c(22000,-100), axes=F, mgp=c(2,0,0), 
			xlab=paste("char -- site", sitenum), ylab="age", pch=16, cex=0.8, type="o")
		axis(1); axis(2); axis(3)

		par(new=T, fig=c(1/3,2/3,2/3,1), mar=c(5.1, 4.1, 2.1, 2.1))
		hist(sitedata$char, xlab="char", mgp=c(2,1,0), main=NULL)

		par(new=T, fig=c(1/3,2/3,1/3,2/3), mar=c(7.1, 4.1, 0.1, 2.1))
		hist(char3, xlab="transformed char", mgp=c(2,1,0), main=NULL)

		par(new=T, fig=c(1/3,2/3,0,1/2), mar=c(4.1, 4.1, 5.1, 2.1))
		plot(liksave ~ lasave, xlim=c(-2,2), xlab=paste("Lambda, opt=",lafit), mgp=c(2,1,0), ylab="Likelihood", type="l")

		par(new=T, fig=c(2/3,1,0,1), mar=c(4.1, 4.1, 2.1, 2.1))
		plot(sitedata$age ~ char3, ylim=c(22000,-100), axes=FALSE, mgp=c(2,0,0), 
			xlab=" transformed char", ylab="age", pch=16, cex=0.8, type="o")
		axis(1); axis(2); axis(3)

	dev.off()

	}
		
}

	

	
	