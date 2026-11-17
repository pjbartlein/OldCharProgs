# presample-bin.R

# presamples, or bins the transformed data into evenly spaced bins with no interpolation  

# paths for input and output .csv files -- modify as appropriate
datapath <- "e:/Projects/GPWG/GPWGv3/data/v3d/"
regionpath <- "e:/Projects/GPWG/GPWGv3/data/v3d/v3h2_regions/"
sitelistfile <- "v3h2_nsa.csv"
queryname <- "v3d"

## set basename and bin structure
basename <- "zt21k"
targbeg <- -60
targend <- 24000
targstep <- 20

# no changes below here

# various path and filenames
sitelistfile <- paste(regionpath,sitelistfile, sep="")
sitelistfile
transcsvpath <- paste(datapath,queryname,"_trans_csv_R/",sep="")
transcsvpath
presampcsvpath <- paste(datapath,queryname,"_presamp_csv_R/",sep="")
presampcsvpath
# if output folder does not exist, create it
dir.create(file.path(datapath, paste(queryname,"_presamp_csv_R/",sep="")), showWarnings=FALSE)

# bin center (target points) definition
targage <- seq(targbeg, targend, by=targstep)

# read list of sites
sitelist <- read.csv(sitelistfile, stringsAsFactors=FALSE)
nsites <- length(sitelist[,1])
print(nsites)

# main loop
for (j in seq(1,nsites)) {
# j=1 # for debugging
  sitenum <- sitelist[j,1]
	sitename <- as.character(sitelist[j,2])
	siteidchar <- as.character(sitenum)
	if (sitenum >= 1) siteid <- paste("000", siteidchar, sep="")
	if (sitenum >= 10) siteid <- paste("00", siteidchar, sep="")
	if (sitenum >= 100) siteid <- paste("0", siteidchar, sep="")
	if (sitenum >= 1000) siteid <- paste(    siteidchar, sep="")

	inputfile <- paste(transcsvpath, siteid, "_trans_influx_",basename,".csv", sep="")

	print(j); print(sitenum); print(sitename); print(inputfile)

	# read the input data
  sitedata <- read.csv(inputfile)
  
  # check to see if all zt values are NA, if so, skip this site
  nonmiss <- na.omit(sitedata$zt)
  length(nonmiss)
  
  if (length(nonmiss) > 0) {
  
    # add a column of 1's for counting
    sitedata$one <- rep(1,length(sitedata[,1]))
    
    #head(cbind(sitedata$est_age, sitedata$zt))
  
    # this definition of bin number seems to match that implicit in presample.f90
    binnum <- as.integer(ceiling((sitedata$est_age-targbeg-(targstep/2))/targstep))+1
    # head(cbind(sitedata$est_age,sitedata$zt,binnum,targage[binnum]), 20)
    
    # get average zt values (and average ages) for the data in each bin
    binave <- tapply(sitedata$zt, binnum, mean)
    binaveage <- tapply(sitedata$est_age, binnum, mean)
    bincount <- tapply(sitedata$one, binnum, sum)
    
    # bin numbers of each bin that had an average (or a single) value
    binsub <- as.numeric(unlist(dimnames(binave)))  
    #head(cbind(targage[binsub],binave,binaveage,bincount),10)
    
    # write output
    presampout <- data.frame(cbind(targage[binsub],binave,bincount))
    presampout <- na.omit(presampout)
    colnames(presampout) <- c("age", "zt", "np")
 
    outputfile <- paste(presampcsvpath, siteid, "_presamp_influx_",basename,"_bw",
      as.character(targstep),".csv", sep="")
    write.table(presampout, outputfile, col.names=TRUE, row.names=FALSE, sep=",")
  }
  
}


