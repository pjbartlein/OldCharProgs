impd <- read.csv("e:\\Projects\\IMPD\\2010\\data\\prop_scar_ann_w104.csv")

hist(impd$PropScar, breaks=20)
plot(impd$Year, impd$PropScar)

# ZINB models
# library(AED)
data(ParasiteCod)
summary(ParasiteCod)

ParasiteCod$fArea <- factor(ParasiteCod$Area)
ParasiteCod$fYear <- factor(ParasiteCod$Year)
I1 <- is.na(ParasiteCod$Intensity) | is.na(ParasiteCod$fArea) |
        is.na(ParasiteCod$fYear) | is.na(ParasiteCod$Length)
ParasiteCod2 <- ParasiteCod[!I1, ]

plot(table(ParasiteCod2$Intensity), ylab="frequencies",
    xlab="observed Intensity values")

# library(pscl)
# library(lmtest)
f1 <- formula(Intensity ~ fArea*fYear + Length | fArea*fYear + Length)
Zip1 <- zeroinfl(f1, dist="poisson", link = "logit", data=ParasiteCod2)
summary(Zip1)

Nb1 <- zeroinfl(f1, dist="negbin", link = "logit", data=ParasiteCod2)
summary(Nb1)

lrtest(Zip1, Nb1)