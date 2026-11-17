cygtrans <- read.csv(file.choose())
summary(cygtrans)
cygtrans2 <- cbind(cygtrans[4],cygtrans[6:16])
plot(cygtrans2, pch=16, cex=0.2)

plot(cygtrans2[1:6], pch=16, cex=0.7)
plot(cygtrans2[c(1:2,7:9)], pch=16, cex=0.7)
plot(cygtrans2[c(1,4:6,10:12)], pch=16, cex=0.7)
plot(cygtrans2[c(1,7:12)], pch=16, cex=0.7)


