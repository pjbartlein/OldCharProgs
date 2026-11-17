insol <- read.csv(file.choose())
# e:\Projects\solar\data\B&L90\file1_90_1Ma.csv

attach(insol)
hist(ECC)
quantile(ECC)
hist(OBL)