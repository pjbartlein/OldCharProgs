# trans_plot.R

# 1-parameter Box-Cox transformation of charcoal quanities for a single site
# (alpha (shift parameter) is specified, lambda (power transformation parameter) is estimated)

# input .csv files should contain at least the variables "est_age" and "quant",
# 	identified by labels in a header row  

# paths for input and output .csv files, and for plots -- modify as appropriate
sitelistcsv <- "e:\\PROJECTS\\GPWG\\data\\files\\v2n\\v2n_sites.csv"
sitecsvpath <- "e:\\Projects\\GPWG\\data\\csv_files\\v2n\\v2n_site\\"
pdfpath <- "e:\\Projects\\GPWG\\data\\plots\\trans_plots\\"

# read list of sites
sitelist <- read.csv(sitelistcsv, stringsAsFactors=FALSE)
nsites <- length(sitelist$ID_SI)

for (j in seq(1,nsites)) {
#for (j in seq(1,1)) {

	sitename <- as.character(sitelist$Sitename[j])
	sitenum <- sitelist$ID_SI[j]
	#print(sitenum)
	siteidchar <- as.character(sitenum)
	if (sitenum >= 1) siteid <- paste("000", siteidchar, sep="")
	if (sitenum >= 10) siteid <- paste("00", siteidchar, sep="")
	if (sitenum >= 100) siteid <- paste("0", siteidchar, sep="")
	if (sitenum >= 1000) siteid <- paste(    siteidchar, sep="")

	inputfile <- paste(sitecsvpath, siteid, "_data.csv", sep="")
	pdfname <- paste(pdfpath, siteid, "_trans.pdf", sep="") 

	print(sitenum)	
	print(sitename)
	print(inputfile)
	print(pdfname)

	# read the input data
        sitedata <- read.csv(inputfile)

        # initial minimax rescaling of data, and addition of "shift" parameter (alpha)
        alpha <- 0.01  # value used in DP1
        #alpha <- 0.5*min(sitedata$quant[sitedata$quant != 0])  # alternative alpha: 0.5 time the smallest nonzero value of quant

        quant2 <- (sitedata$influx2-min(sitedata$influx2))/(max(sitedata$influx2)-min(sitedata$influx2))+alpha
        ###quant2 <- (sitedata$conc2-min(sitedata$conc2))/(max(sitedata$conc2)-min(sitedata$conc2))+alpha

        # maximum likelihood estimation of lambda
        # derived from the boxcox.R function in the Venables and Ripley MASS library included in R 2.6.1

        npts <- 201 # number of estimates of lambda
        y <- quant2
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
        cbind(sitename,liksave[which.max(liksave)],lasave[which.max(liksave)])
        laopt=as.character(lasave[which.max(liksave)])
        lafit <- lasave[which.max(liksave)]

        # Box-Cox transformation of data
        if (laopt == 0.0) tquant <- log(quant2) else tquant <- (quant2^lafit - 1)/lafit

        # if comparing multiple records, do minimax rescaling
        tquant <- (tquant - min(tquant))/(max(tquant)-min(tquant))

        # output data and plots (if file errors occur, .csv and .pdf files are probably already open)

        # plot the original data, histograms, and transformed data -- use one of the following and see dev.off() below
        pdf(file=pdfname,width=8.5, height=11.0)  # create a .pdf file of plots
        # windows(width=8.5, height=11.0)  # plot to the windows device (within R)

        # plot raw data
        plotx <- sitedata$est_age[sitedata$est_age != -9999]
        ploty <- quant2[sitedata$est_age != -9999]
        par(fig=c(0,1,.70,.95), mar=c(3.1, 4.1, 2.1, 4.1))
        #plot(sitedata$influx2 ~ sitedata$est_age, xlim=c(20000,-100), axes=F, mgp=c(2,0,0),
        plot(ploty ~ plotx, xlim=c(max(sitedata$est_age),-100), axes=F, mgp=c(2,0,0),
                main=paste("Box-Cox Transformation -- ", sitename, "Site ", siteidchar), font.main=1, lab=c(10,5,5),
                ylab="influx2", xlab="Age", cex.lab=0.8, pch=16, cex=0.5, type="o")
        ####plot(ploty ~ plotx, xlim=c(max(sitedata$est_age),-100), axes=F, mgp=c(2,0,0),
		####                main=paste("Box-Cox Transformation -- ", sitename, "Site ", siteidchar), font.main=1, lab=c(10,5,5),
        ####        ylab="conc2", xlab="Age", cex.lab=0.8, pch=16, cex=0.5, type="o")
        
        #axis(1, cex.axis=0.8, xaxp=c(0,22000,22)); axis(2, cex.axis=0.8); axis(4, cex.axis=0.8)
	    axis(1, cex.axis=0.8); axis(2, cex.axis=0.8); axis(4, cex.axis=0.8)


        # plot transformed data
        plotx <- sitedata$est_age[sitedata$est_age != -9999]
        ploty <- tquant[sitedata$est_age != -9999]
        par(new=T, fig=c(0,1,.4,.65), mar=c(3.1, 4.1, 2.1, 4.1))
	    #plot(tquant ~ sitedata$est_age, xlim=c(20000,-100), axes=F, mgp=c(2,0,0), lab=c(10,5,5),
        plot(ploty ~ plotx,xlim=c(max(sitedata$est_age),-100), axes=F, mgp=c(2,0,0), lab=c(10,5,5),
                ylab="transformed influx2", xlab="Age", cex.lab=0.8, pch=16, cex=0.5, type="o")
        ####plot(ploty ~ plotx,xlim=c(max(sitedata$est_age),-100), axes=F, mgp=c(2,0,0), lab=c(10,5,5),
        ####        ylab="transformed conc2", xlab="Age", cex.lab=0.8, pch=16, cex=0.5, type="o")
        #axis(1, cex.axis=0.8, xaxp=c(0,22000,22)); axis(2, cex.axis=0.8); axis(4, cex.axis=0.8)
        axis(1, cex.axis=0.8); axis(2, cex.axis=0.8); axis(4, cex.axis=0.8)

        # plot histograms of raw and transformed data, and likelhood profile
        par(new=T, fig=c(0,.4,.0,.35), mar=c(7.1, 4.1, 2.1, 2.1))
        hist(sitedata$influx2, xlab="influx2", cex.lab=0.8, cex.axis=0.8, mgp=c(2,1,0), main=NULL)
        ####hist(sitedata$conc2, xlab="conc2", cex.lab=0.8, cex.axis=0.8, mgp=c(2,1,0), main=NULL)

        par(new=T, fig=c(.3,.7,.0,.35), mar=c(7.1, 4.1, 1.1, 2.1))
        hist(tquant, xlab="trans influx2", cex.lab=0.8, cex.axis=0.8, mgp=c(2,1,0), main=NULL)
        ####hist(tquant, xlab="trans conc2", cex.lab=0.8, cex.axis=0.8, mgp=c(2,1,0), main=NULL)

        par(new=T, fig=c(.6,1.0,.0,.35), mar=c(7.1, 4.1, 1.1, 2.1))
        plot(liksave ~ lasave, xlim=c(-2,2), xlab=paste("Lambda, opt=",format(lafit,digits=3)),
                cex.lab=0.8, cex.axis=0.8, mgp=c(2,1,0), ylab="Likelihood", type="l")
                title(xlab=paste("Alpha =", format(alpha,digits=3)), cex.lab=0.8)

        dev.off() # turn off the .pdf device (comment this line out if plotting to the windows device


}	