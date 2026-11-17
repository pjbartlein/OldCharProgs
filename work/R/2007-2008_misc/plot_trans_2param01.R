
csvpath <- "E:\\Projects\\igbpfire\\charsum\\data\\current_hs2\\sites\\site_csv\\"
pdfpath <- "E:\\Projects\\igbpfire\\charsum\\data\\current_hs2\\sites\\transplot\\"

alpha <- c(0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1.0,2.0,5.0,
			10.0,20.0,50.0,100.0,200.0,500.0,1000.0)
nalpha <- 22


#alpha=.01
#nalpha=1
la <- seq(-2, 2, by=0.04)
nlam <- length(la)


j=221

#for (j in seq(1,700)) {	
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
		sitedata2 <- read.csv(csvfile) 
		#sitedata2$quant <- sitedata2$quant
		#print(csvfile) 
	
		lasave <- matrix(nrow=nlam, ncol=nalpha)
		alphasave <- matrix(nrow=nlam, ncol=nalpha)
		liksave <- matrix(nrow=nlam, ncol=nalpha)
		
		for (ii in 1:nalpha) {
			for (i in 1:nlam) {
				
				quant2 <- sitedata2$quant #+alpha[ii]
				#quant2 <- (sitedata2$quant-min(sitedata2$quant))/(max(sitedata2$quant)-min(sitedata2$quant))+alpha[ii]
				y <- quant2
				n <- length(y)
				logy <- log(y)
				ydot <- exp(mean(log(y+alpha[ii])))
								
				#if (la[i] != 0.0) yt <- ((y^la[i])-1)/la[i] else yt <- logy*(1+(la[i]*logy)/2*(1+(la[i]*logy)/3*(1+(la[i]*logy)/4))) 
				if (la[i] != 0.0) yt <- (((y+alpha[ii])^la[i])-1)/la[i] else yt <- log(y+alpha[ii])
				
				#zt <- yt/(ydot^(la[i]-1)) # one parameter		
				zt <- (((y+alpha[ii])^la[i])-1)/(la[i]*(ydot^(la[i]-1)))
				#jac <- ((la[i]-1)*sum(log(y)))
				#zt <- yt/jac^(1/n)
				loglik <- -n/2*log(sum((zt - mean(zt))^2 ))
				#loglik <- -n/2*log((sum((zt - mean(zt))^2)/n))
				#loglik <- -n/2*log(sum((yt-mean(yt))^2))
				#loglik <- -n/2*log(sum((yt - mean(yt))^2 ))+((la[i]-1)*sum(log(y)))
				
				#xmat <- rep(1, n) 
				#xmat <- as.matrix(xmat)
				#beta <- solve(crossprod(xmat), crossprod(xmat, yt))
  				#ss <- sum((drop(yt) - drop(xmat %*% beta))^2) 
				#loglik <- (n/2) * log(ss) - ((la[i] - 1) * sum(log(y))) 
				
				lasave[i,ii] <- la[i]
				alphasave[i,ii] <- alpha[ii]
				liksave[i,ii] <- loglik
				#print(cbind(i,ii,lasave[i,ii],alphasave[i,ii],liksave[i,ii]))
				
				
				#ss <- sum((drop(yt) - drop(xmat %*% beta))^2)
				#  if(lik.method == "ML")
    			# neglik <- (n/2) * log(ss) - ((lambda - 1) * sum(log(data))) 
			} 
		}
		#print(cbind(lasave,liksave,alphasave))
		cbind(sitenum,liksave[which.max(liksave)],lasave[which.max(liksave)],alphasave[which.max(liksave)])
		alphaopt=as.character(alphasave[which.max(liksave)])
		alphafit <- alphasave[which.max(liksave)]
		laopt=as.character(lasave[which.max(liksave)])
		lafit <- lasave[which.max(liksave)]
		quant2 <- sitedata2$quant+alphafit
		if (laopt == 0.0) quant3 <- log(quant2) else quant3 <- (quant2^lafit - 1)/lafit


		#pdf(file=pdfname ,width=11.0, height=8.5)
		windows(width=11.0, height=8.5)
		
		par(fig=c(0,1/3,0,1), mar=c(4.1, 4.1, 2.1, 2.1))
		plot(sitedata2$est_age ~ sitedata2$quant, ylim=c(22000,-100), axes=F, mgp=c(2,0,0), 
			xlab=paste("quant -- site", sitenum), ylab="age", pch=16, cex=0.8, type="o")
		axis(1); axis(2); axis(3)

		par(new=T, fig=c(1/3,2/3,2/3,1), mar=c(5.1, 4.1, 2.1, 2.1))
		hist(sitedata2$quant, xlab="quant", mgp=c(2,1,0), main=NULL)

		par(new=T, fig=c(1/3,2/3,1/3,2/3), mar=c(7.1, 4.1, 0.1, 2.1))
		hist(quant3, xlab="transformed quant", mgp=c(2,1,0), main=NULL)

		par(new=T, fig=c(1/3,2/3,0,1/2), mar=c(4.1, 4.1, 5.1, 2.1))
		
		image(la,log(alpha),-1.0*liksave, xlab=paste("Lambda, opt=",lafit), 
			ylab=paste("log(Alpha), opt=",alphafit), mgp=c(2,1,0), col=rainbow(50))
		#contour(la,log(alpha),liksave, xlab=paste("Lambda, opt=",lafit), 
		#	ylab=paste("log(Alpha), opt=",alphafit), mgp=c(2,1,0), nlevels=40, add=TRUE)

		
		par(new=T, fig=c(2/3,1,0,1), mar=c(4.1, 4.1, 2.1, 2.1))
		plot(sitedata2$est_age ~ quant3, ylim=c(22000,-100), axes=FALSE, mgp=c(2,0,0), 
			xlab=" transformed quant", ylab="age", pch=16, cex=0.8, type="o")
		axis(1); axis(2); axis(3)

	#dev.off()

	}
		
#}

	

	
	