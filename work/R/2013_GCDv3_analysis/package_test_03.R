# compare paleofire with Fortran and locfit versions

## Install from Github
#install.packages("devtools")
library(devtools)
install_github("paleofire","paleofire")
library(paleofire)


# select sites -- ENA
enaID <- pfSiteSel(Latlim=c(30,70),Longlim=c(-100,-50))
plot(enaID)
length(enaID$SiteNames)

# summary of site list
ena_list <- summary(enaID)
head(ena_list)

# write it out to compare with ena_sites_a-ld.csv
sitecsv <- "package_test_02_ena_sitelist.csv"
write.table(ena_list, sitecsv, row.names=FALSE, sep=",")

# transform data
ena_trans <- pfTransform(enaID, method=c("MinMax","Box-Cox","MinMax","Z-Score"),
  BasePeriod=c(200,4000),QuantType="INFL")

# downcore plots with transformed series (pdf file:
pfDiagnostic(enaID, method=c("Box-Cox","Z-Score"),BasePeriod=c(200,4000))

## Compositing: Using the locfit package equivalent procedure to Daniau et al. 2012
ena_fit <- pfCompositeLF(ena_trans, tarAge=seq(-60,8000,20), binhw=10, hw=200, nboot=1000)
plot(ena_fit)

## And save
csvfile <- "v3h_ena_sites_paleofire_zt4k_bw20_200_1000.csv"
names(ena_fit$Result) <- c("AGE", "pfLocFit",  "pfMEAN(of_boot)", "pfCI05", "pfCI95" )
write.table(ena_fit$Result, csvfile, row.names=FALSE, sep=",")


# select sites -- globe
v3ID <- pfSiteSel(Latlim=c(-90,90),Longlim=c(-180,180))
plot(v3ID)
length(v3ID$SiteNames)

# get rid of site 882
v3ID$SiteNames <- v3ID$SiteNames[v3ID$SitesIDS != 882]
v3ID$SitesIDS <- v3ID$SitesIDS[v3ID$SitesIDS != 882]
plot(v3ID)
length(v3ID$SiteNames)

# summary of site list
v3_list <- summary(v3ID)
head(v3_list)

# write it out to compare with ena_sites_a-ld.csv
sitecsv <- "package_test_02_v3_sitelist.csv"
write.table(v3_list, sitecsv, row.names=FALSE, sep=",")

# transform data
v3_trans <- pfTransform(v3ID, method=c("MinMax","Box-Cox","MinMax","Z-Score"),
  BasePeriod=c(200,4000),QuantType="INFL")

# downcore plots with transformed series (pdf file:
#pfDiagnostic(v3ID, method=c("Box-Cox","Z-Score"),BasePeriod=c(200,4000))

## Compositing: Using the locfit package equivalent procedure to Daniau et al. 2012
v3_fit <- pfCompositeLF(v3_trans, tarAge=seq(-60,8000,20), binhw=10, hw=200, nboot=1000)
plot(v3_fit)

## And save
csvfile <- "v3h_v3_sites_paleofire_zt4k_bw20_200_1000.csv"
names(v3_fit$Result) <- c("AGE", "pfLocFit",  "pfMEAN(of_boot)", "pfCI05", "pfCI95" )
write.table(v3_fit$Result, csvfile, row.names=FALSE, sep=",")

