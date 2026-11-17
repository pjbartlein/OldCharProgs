# Extract relevant data from GCD v.1 Clim. Dyn. paper
# This script is based on Simon's QUAVIDA\getData.r script
# Should adjust this to use site nos ...


require(RODBC)
ch <- odbcConnectAccess("charcoal-06-Jan-09.mdb")
#ch <- odbcConnectAccess("GCDv1_n405jrm.mdb")
sites <- sqlFetch(ch, "SITE")
samples <- sqlFetch(ch, "SAMPLE")
charcoal <- sqlFetch(ch, "CHARCOAL")
odbcClose(ch)

ents <- (unique(sites$ID_SITE))
snames <- as.character(unique(sites$SITE_NAME))

csvpath <- "c:\\data\\projects\\R\\ppp\\out\\"
for (j in 1:length(ents)) {
        sitenum <- j
      	sina <- snames[j]
        siteidchar <- as.character(sitenum)
      	if (sitenum >= 1) siteid <- paste("000", siteidchar, "_", sep="")
      	if (sitenum >= 10) siteid <- paste("00", siteidchar, "_", sep="")
      	if (sitenum >= 100) siteid <- paste("0", siteidchar, "_", sep="")
      	print(siteid)
        csvfile <- paste(csvpath, siteid, sina, ".csv", sep="")
      	samps <- samples[which(samples$ID_SITE==ents[j]),]
        #charc <- charcoal[which(charcoal$ID_SITE==ents[j]),]
        write.table(samps, file = csvfile, sep=",", row.names=FALSE)
        }
odbcClose(ch)
####################################

