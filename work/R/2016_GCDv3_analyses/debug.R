debug <- file("debug.txt", "w")

writeLines(siteidchar, con = debug, sep = "\n")
writeLines(inputfile, con = debug, sep = "\n")