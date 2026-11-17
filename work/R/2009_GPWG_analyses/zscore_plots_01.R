# zscore_plots.R

# plots of z-scores for various base periods

# paths for input and output .csv files, and for plots -- modify as appropriate
sitelistcsv <- "\\PROJECTS\\GPWG\\data\\files\\v2l\\v2l_sites_combined.csv"
sitecsvpath <- "\\Projects\\GPWG\\data\\csv_files\\v2l_trans_influx2\\"
pdfpath <- "\\Projects\\GPWG\\data\\plots\\zscore_plots\\"

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

        inputfile <- paste(sitecsvpath, siteid, "_trans_influx2.csv", sep="")
        pdfname <- paste(pdfpath, siteid, "_zscores.pdf", sep="")

        print(sitenum)
        print(sitename)
        print(inputfile)
        print(pdfname)

        # read the input data
        sitedata <- read.csv(inputfile)

        # plot the zscores -- use one of the following and see dev.off() below
        pdf(file=pdfname,width=8.5, height=11.0)  # create a .pdf file of plots
        # windows(width=8.5, height=11.0)  # plot to the windows device (within R)

        # plot raw data
        plotmax <- max(c(sitedata$ztall, sitedata$zt21k, sitedata$ztpi4, sitedata$ztpi2))
        plotmin <- min(c(sitedata$ztall[sitedata$ztall != -9999.], sitedata$zt21k[sitedata$zt21k != -9999.],
                        sitedata$ztpi4[sitedata$ztpi4 != -9999.], sitedata$ztpi2[sitedata$ztpi2 != -9999.]))
        par(fig=c(0,1,.70,.95), mar=c(3.1, 4.1, 2.1, 4.1))
        plot(sitedata$ztall ~ sitedata$est_age, xlim=c(max(sitedata$est_age),-100),
                        ylim=c(plotmin,plotmax), axes=F, mgp=c(2,0,0),
                main=paste("Z-Scores -- ", sitename, "Site ", siteidchar), font.main=1, lab=c(10,5,5),
                ylab="Z-Score", xlab="Age", cex.lab=0.8, pch=16, cex=0.5, type="o")
                axis(1, cex.axis=0.8); axis(2, cex.axis=0.8); axis(4, cex.axis=0.8)
                lines(sitedata$est_age, sitedata$ztpi2, pch=16, cex=0.5, type="o", col="gray")
                lines(sitedata$est_age, sitedata$zt21k, pch=16, cex=0.5, type="o", col="purple")
                lines(sitedata$est_age, sitedata$ztpi4, pch=16, cex=0.5, type="o", col="blue")
                abline(0.0, 0.0)

        dev.off() # turn off the .pdf device (comment this line out if plotting to the windows device


}

