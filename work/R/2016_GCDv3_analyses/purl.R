library(knitr)

filename <- "trans-and-zscore_06"
purl(paste(filename, ".Rmd", sep=""), paste(filename, "_code.R", sep=""), documentation=0)

