library(mzR)
msfile<-'000_select.mzXML'
my.data<-openMSfile(msfile)

## figure 11.1

plist<-peaks(my.data,7)

tiff("figure11_1.tif",width=5,height=5,units="in",res=1200)
plot(plist, type="h", xlab="m/z",ylab="Intensity",main="Scan 7")
dev.off()

pl<-peaks(my.data,7)

## figure 11.3

topnum<-30
pl.selected<-pl[pl[,2]%in% head(sort(pl[,2],decreasing=T),topnum),]

tiff("figure11_3.tif",width=5,height=5,units="in",res=1200)
par(mfrow=c(2,1),mai=c(0.55,0.2,0.1,0.1))
plot(pl,type="h",xlim = c(300,1000)) 
mtext("A", side = 3, line = -2, adj = 0.01, cex = 1.3)
plot(pl.selected,type="h",xlim = c(300,1000)) 
mtext("B", side = 3, line = -2, adj = 0.01, cex = 1.3)
dev.off()
