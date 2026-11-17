install.packages("devtools")
library(devtools)

install_github("paleofiRe","oblarquez")
library(paleofiRe)

## Simple test

ID=pfSiteSel(Region="ASIA")
plot(ID)


# Select sites in Eastern North America:
ID=pfSiteSel(Region="ENA0")

# Name of selected sites
ID$SiteNames

# Site extraction 
A=pfExtract(ID)

# Site 1 Serie Example plot
plot(A[A[,1]==ID$SitesIDS[1],3],A[A[,1]==ID$SitesIDS[1],4],type="l",main=ID$SiteNames[1])

# Plot the map of selected sites
pfPlotMap(ID)
fifdev.off()

plot(ID, type="chronology"


ls(2)