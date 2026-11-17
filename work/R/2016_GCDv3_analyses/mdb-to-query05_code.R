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

# site .csv file
csvpath <- "/Projects/GPWG/GPWGv3/GCDv3Data/v3i/v3i_query/"
sitefile <- "v3i_SITES.csv"
write.csv(site_query, paste(csvpath, sitefile, sep=""), row.names=FALSE)

# data .csv file
csvpath <- "/Projects/GPWG/GPWGv3/GCDv3Data/v3i/v3i_query/"
datafile <- "v3i_DATA.csv"
write.csv(data_query, paste(csvpath, datafile, sep=""), row.names=FALSE)
