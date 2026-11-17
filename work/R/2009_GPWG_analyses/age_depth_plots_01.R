# age_depth_plots.R

# paths for input and output .csv files, and for plots -- modify as appropriate
sitelistcsv <- "\\PROJECTS\\GPWG\\data\\files\\v2n\\v2n_sites.csv"
sitecsvpath <- "\\Projects\\GPWG\\data\\csv_files\\v2n\\v2n_site\\"
pdfpath <- "\\Projects\\GPWG\\data\\plots\\age_depth_plots\\"

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
	pdfname <- paste(pdfpath, siteid, "_age_depth_plot.pdf", sep="")
	print(sitenum)	
	print(sitename)
	print(inputfile)
	print(pdfname)

	# read the input data
      sitedata <- read.csv(inputfile)       

      # output data and plots (if file errors occur, .csv and .pdf files are probably already open)

      # plot the original data, histograms, and transformed data -- use one of the following and see dev.off() below
      pdf(file=pdfname) #,width=8.5, height=11.0)  # create a .pdf file of plots
      #windows(width=8.5, height=11.0)  # plot to the windows device (within R)

      # plot raw data
      plot(sitedata$depth ~ sitedata$est_age, 
      	main=paste(sitename, "  ", siteidchar), font.main=1, lab=c(10,5,5),
            ylab="Depth", xlab="Age", cex.lab=0.8, pch=16, cex=0.5)
        

      dev.off() # turn off the .pdf device (comment this line out if plotting to the windows device


}	

