# mdb-to-query.R
# directly reads an Access database and creates site and data query files

# load RODBC library and connect to the database
library(RODBC)
dbname <- "GCDv03_Marlon_et_al_2015.mdb"
gcdv3.db <- odbcConnect(dbname)
odbcGetInfo(gcdv3.db)

# check for existence of database site and data views
sqlTables(gcdv3.db, tableName="ALL_BART_SITES", tableType="VIEW")
sqlColumns(gcdv3.db, "ALL_BART_SITES")$COLUMN_NAME
sqlTables(gcdv3.db, tableName="ALL_BART_DATA", tableType="VIEW")
sqlColumns(gcdv3.db, "ALL_BART_DATA")$COLUMN_NAME

# site query
site_query <- sqlFetch(gcdv3.db, "ALL_BART_SITES")
names(site_query)[3] <- "LATITUDE"; names(site_query)[4] <- "LONGITUDE"
site_query$SITE_NAME <- as.character(site_query$SITE_NAME)
head(site_query)
str(site_query)

#data query
data_query <- sqlFetch(gcdv3.db, "ALL_BART_DATA")
names(data_query)[4] <- "DEPTH"
head(data_query)
str(data_query)

# close the database
odbcClose(gcdv3.db)

# loop over sites

miss <- -9999.0
j <- 1

  nsamp <- 0
  sitedata <- data_query[data_query$ID_SITE == j, ]
  nsamp <- length(sitedata$ID_SITE)
  head(sitedata)
  str(sitedata)
  
  if (nsamp > 0) {
    jchar <- as.character(j)
    nsampchar <- as.character(nsamp)
  
    # local variables
    depth <- sitedata$DEPTH; age <- sitedata$EST_AGE; quant <- sitedata$QUANTITY
    depth[is.na(depth)] <- miss
    age[is.na(age)] <- miss
    quant[is.na(quant)] <- miss
    
    thickness <- rep(miss, nsamp); dep_time <- rep(miss, nsamp); sed_rate <- rep(miss, nsamp)
    unit_dep_time <- rep(miss, nsamp)
    xst_level <- as.character(sitedata[1,9])
    
    # sed rate and deposition time
    
    # first (top) sample
    if (depth[1] != miss && depth[2] != miss) {
      thickness[1] <- (depth[2] - depth[1])*100.0 # meters to cm (depth in m, influx and conc in cm)
      dep_time[1] <- age[2] - age[1]
      if (dep_time[1] > 0.0) sed_rate[1] <- thickness[1]/dep_time[1]
      if (sed_rate[1] != miss) unit_dep_time[1] <- 1.0/sed_rate[1]
    }
    # samples 2 to nsamp-1
    for (i in 2:(nsamp-1)) {
      if (depth[1] != miss && depth[2] != miss) {
        thickness[i] <- (depth[i+1] - depth[i])*100.0 
        dep_time[i] <- ((age[i+1] + age[i])/2.0) - ((age[i] + age[i-1])/2.0)
        if (dep_time[i] > 0.0) sed_rate[i] <- thickness[i]/dep_time[i]
        if (sed_rate[i] != miss) unit_dep_time[i] <- 1.0/sed_rate[i] 
      }
    }
    # last (bottom) sample
    if (depth[nsamp-1] != miss  && depth[nsamp] != miss) {
      thickness[nsamp] <- thickness[nsamp-1] # replicate thickness
      dep_time[nsamp] <- age[nsamp] - age[nsamp-1]
      sed_rate[nsamp] <- sed_rate[nsamp-1] # replicate sed_rate
      unit_dep_time[nsamp] <- unit_dep_time[nsamp-1]
    }
    
    # counts of missing values
    depth_count <- 0; age_count <- 0; quant_count <- 0; sed_rate_count <- 0; sed_rate_flag <- 1
    depth_count <- sum(depth != miss)
    age_count <- sum(age != miss)
    quant_count <- sum(quant != miss)
    sed_rate_count <- sum(sed_rate != miss)
    if (sed_rate_count != nsamp) sed_rateflag = 0
    
    # check for age or depth reversals, and zero or negative sed rates (in nonmissing data)
    depth_reversal <- 0; age_reversal <- 0; sed_rate_zeroneg <- 0         
    for (i in 2:nsamp) {
      if (age[i] != miss && age[i-1] != miss && age[i] <= age[i-1]) age_reversal=1
      if (depth[i] != miss && depth[i-1] != miss) {
        if (depth[i] <= depth[i-1]) depth_reversal=1
      } 
    }
    for (i in 2:nsamp) {
      if (sed_rate[i] != miss && sed_rate[i] <= 0.0) sed_rate_zeroneg=1
    }
    
    # alternative quantities
    
    conc <- rep(miss, nsamp); influx <- rep(miss, nsamp)
    influx_source <- rep("none", nsamp) ; conc_source <- rep("none", nsamp)
    
    # select case based on 1st_level
    
    if (xst_level == "INFL")          # adopt influx values as they are, calculate concentration
    {  
      influx <- quant
      influx_source <- "data"
      if (influx != miss && unit_dep_time != miss && sed_rate != 0.0) {
        conc <- influx * unit_dep_time
        conc_source <- "calculated from influx "
      } else {
        conc <- quant
        conc_source <- "copied from quant "
      }
    }
    
    else if (xst_level == "CONC")     # calculate influx, adopt conc values as they are
    {
      conc <- quant
      conc_source <- "data"
      if (conc != miss && sed_rate != miss && sed_rate != 0.0) {
        influx <- quant * sed_rate
        influx_source <- "calculated from conc "
      } else {
        influx <- quant
        influx_source <- "copied from quant "
      }  
    }
    
    else if (xst_level == "C0P0")     # assume quantity is concentration like
    {
      conc <- quant
      conc_source <- "C0P0"
      if (sed_rate != miss && sed_rate != 0.0) {
        influx <- quant * sed_rate
        influx_source <- "calculated from C0P0 (conc) "
      } else {
        influx <- quant
        influx_source <- "copied from quant "
      }    
    }
    
    else if (xst_level == "SOIL")     # just copy
    {
      conc <- quant
      conc_source <- "copied from quant "
      influx <- quant
      influx_source <- "copied from quant "
    }
    else if (xst_level == "OTHE")     # just copy
    {
      conc <- quant
      conc_source <- "copied from quant "
      influx <- quant
      influx_source <- "copied from quant "
    }
    else 
    {
      conc <- quant
      conc_source <- "copied from quant "
      influx <- quant
      influx_source <- "copied from quant "
    }
  }
  
  # check for influx == 0.0 everywhere
  nzero <- 0
  nzero <- sum(influx != 0.0)
  if (nzero == 0) {
    writeLines(paste("**** has no non-zero influx values", sep=" "), con = debugfile, sep = "\n")
  }
  
  # .csv out
  if (nsamp > 0 && nzero > 0) {
    
    # get siteid string
    siteidchar <- as.character(j)
    if (j >= 1) siteid <- paste("000", siteidchar, sep="")
    if (j >= 10) siteid <- paste("00", siteidchar, sep="")
    if (j >= 100) siteid <- paste("0", siteidchar, sep="")
    if (j >= 1000) siteid <- paste(    siteidchar, sep="")
    sitehdr <- paste("site", siteid, sep="")
    
    # assemble output data and write it out
    samplenum <- seq(1:nsamp)
    outdata <- data.frame(samplenum,sitedata$ID_SAMPLE, depth, age, sed_rate, quant, conc, 
                          influx, xst_level, conc_source, influx_source)
    names(outdata) <- c(sitehdr, "id_sample", "depth", "est_age", "sed_rate", "quant", "conc",
                        "influx", "xst_level", "conc_source", "influx_source" )
  }

save.image("E:/Projects/GPWG/GPWGv3/R/mdb-to-csv_06_run1.RData")
