#Check GCD query output for age reversals and negative sed rates
#"ID_SITE","PREF_UNITS","ID_SAMPLE","Expr1","EST_AGE","QUANTITY","ID_DEPO_CONTEXT","ID_STATUS","1ST_LEVEL"

gcd <- read.csv(file.choose()) # raw charcoal data
head(gcd) # make sure first row contains labels
sitelist <- read.csv(file.choose()) # site list
head(sitelist) # make sure first row contains labels
nrow(gcd) # print the number of rows in the query
sids <- unique(gcd$ID_SITE)  # list of site ids
length(sids) # number of sites in the file

for (j in 1:length(sids)) {
    current.sid <- sids[j]
    #print(paste("Checking DEPTHS for site", current.sid))
    tmpDepths <- gcd[which(gcd$ID_SITE==sids[j]),4]
       if (is.na(tmpDepths)) next else
            tDL <- length(tmpDepths)
            for (i in 1:(tDL-1)) {
                if (tmpDepths[i] >= tmpDepths[i+1])
                  print(paste("Depths are OUT OF ORDER", current.sid))
            }
            rm(tmpDepths)
} 

## Same for Ages 
 for (j in 1:length(sids)) {
    current.sid <- sids[j]
    #print(paste("Checking AGES for site", current.sid))
    tmpAges <- gcd[which(gcd$ID_SITE==sids[j]),5]
       if (is.na(tmpAges)) next else
            tDL <- length(tmpAges)
            for (i in 1:(tDL-1)) {
                if (tmpAges[i] >= tmpAges[i+1])
                  print(paste("Ages are OUT OF ORDER", current.sid))
            }
            rm(tmpAges)
} 

# Check # of sites in data file against # of sites in site list file
write.csv(sids, file="sites.in.data.query.csv")
nrow(sitelist) # number of rows in site list query
length(sids) # of sites in data query                                     
# these should match -- investigate discrepancies
 x<-c(sids,sitelist$ID_SITE)
 table(x)
 # find the table value(s) that is not 2... those sites are only in one file