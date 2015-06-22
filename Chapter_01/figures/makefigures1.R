library(Cairo)
furin.file<- 'furin_data.csv'
my.data<-read.csv(furin.file,sep="\t")

CairoTIFF("figure1_1.tif",width=5,height=5,units="in",res=1200)

#par(oma=c(1,1,1,1))
plot(my.data$Naive.KO.1,ylab="Data",main="Single column from a dataframe",pch=20,cex=0.5)
#plot(my.data)

dev.off()

CairoTIFF("figure1_2.tif",width=5,height=5,units="in",res=1200)

#par(mai=c(0.1,0.1,0.1,0.1),omi=c(0.1,0.1,0.1,0.1))
plot(my.data[3:6],main="Correlation of multiple columns",pch=20,cex=0.2)
#plot(my.data)

dev.off()

CairoTIFF("figure1_3.tif",width=5,height=5,units="in",res=1200)

b<-rnorm(1000)
hist(b)

dev.off()
