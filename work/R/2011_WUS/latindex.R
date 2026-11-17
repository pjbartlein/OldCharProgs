latindex <- read.csv(file.choose())

newlat <- seq(32, 48, by=1)
newseq <- approx(latindex$lat,latindex$index,newlat, rule=2)

outindex <- cbind(newseq$x,newseq$y)
colnames(outindex) <- c("newlat","newseq")
outindex

write.csv(outindex, file.choose(), row.names=FALSE)