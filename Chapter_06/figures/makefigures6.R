#datapath <- '/home/data/work/Courses/Microarray/GSE29797/'
setwd(datapath)
library(affy)
dat<-ReadAffy()

library(Cairo)

CairoTIFF("figure6_1.tif",width=5,height=5,units="in",res=1200)

image(dat[,1])

dev.off()


