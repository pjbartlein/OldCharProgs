# runs knitr and markdown
library(knitr)
library(markdown)

filename <- "regr01"
Rmd.file <- paste(filename, ".Rmd", sep="")
md.file <- paste(filename, ".md", sep="")
html.file <- paste(filename, ".htm", sep="")

# knit the file
knit(Rmd.file, opts_chunk$set(comment=NA))

# convert .md to .html using custom style sheet
markdownToHTML(file=md.file,output=html.file, stylesheet='projector.css')






