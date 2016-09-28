library(mzR)
msfile<-'000_select.mzXML'
my.data<-openMSfile(msfile)

## figure 11.1

plist<-peaks(my.data,7)

tiff("figure11_1.tif",width=5,height=5,units="in",res=1200)
plot(plist, type="h", xlab="m/z",ylab="Intensity")
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

## figure 11.4

topnum<-15
pl.selected<-pl[pl[,2]%in% head(sort(pl[,2],decreasing=T),topnum),]

aa.mass<-read.csv("aa_mass.csv",row.names=1)

peakdiff.aa<-matrix(rep('-',topnum^2),ncol=topnum,nrow=topnum)

topnum<-150
prec<-0.04
pl.top<-pl[pl[,2]%in% head(sort(pl[,2],decreasing=T),topnum),1]
peakdiff<-outer(pl.top,pl.top,'-')
peakdiff.aa<-matrix(rep('-',topnum^2),ncol=topnum,nrow=topnum)
aa.mass<-read.csv("aa_mass.csv",row.names=1)
ystep<-max(pl[,2])/25

tiff("figure11_4.tif",width=5,height=5,units="in",res=1200)
ypos<-ystep*4
plot(pl,type="h",xlab="m/z",ylab="Intensity",main="Scan 7")
for (aa in row.names(aa.mass)){
  peakdiff.aa[abs(peakdiff-aa.mass[aa,])<prec]<-aa              
  aa.match<-which(abs(peakdiff-aa.mass[aa,])<prec,arr.ind=T)
  if(length(aa.match)){
    for (i in 1:dim(aa.match)[1]){
      segments(pl.top[aa.match[i,2]], ypos,pl.top[aa.match[i,1]], ypos,lwd=2)
#      text(pl.top[aa.match[i,2]]+60,ypos+1800000,aa)
    }
  }
  ypos<-ypos+ystep
}
text(850,57035581,"V")
text(970,54035581,"Y")
text(1140,37535581,"K")
dev.off()
    
## figure 11.5

library(MSnbase)
mztab <- 'F063721.dat-mztab.txt'
qnt <- readMzTabData(mztab, what = "PEP")
sampleNames(qnt)<-c("TMT6.126", 
                    "TMT6.127", 
                    "TMT6.128", 
                    "TMT6.129", 
                    "TMT6.130", 
                    "TMT6.131")
qntS <- normalise(qnt, "sum")
qntV <- normalise(qntS, "vsn")
qntV2 <- normalise(qnt, "vsn")


tiff("figure11_5.tif",width=5,height=5,units="in",res=1200)
par(mfrow=c(2,2),mai=c(0.12,0.5,0.5,0.1))
boxplot(exprs(qnt),main="Raw data",names=NA)
mtext("A", side = 3, line = 0, adj = 0.01, cex = 1.1)
boxplot(exprs(qntS),main="Sum normalization",names=NA)
mtext("B", side = 3, line = 0, adj = 0.01, cex = 1.1)
boxplot(exprs(qntV2),main="Variance stabilization\nof raw data",names=NA)
mtext("C", side = 3, line = 0, adj = 0.01, cex = 1.1)
boxplot(exprs(qntV),main="Variance stabilization\nof sum",names=NA)
mtext("D", side = 3, line = 0, adj = 0.01, cex = 1.1)
dev.off()

## figure 11.6

protqnt <- combineFeatures(qnt, groupBy = fData(qnt)$accession, fun = sum)
protexp<-exprs(protqnt)
row.names(protexp)[400:404]<-c("PYGM_RABIT", 
                               "TRYP_PIG", 
                               "ENO1_YEAST", 
                               "ALBU_BOVIN", 
                               "CYC_BOVIN")

tiff("figure11_6.tif",width=5,height=5,units="in",res=1200)
par(mai=c(0.4,0.8,0.4,0.1))
matplot(t(protexp[400:404,]), type = "b", 
        col="black",
        lty=c(1,2,3,4,5),
        ylab="Protein abundance",
        main="Abundance of selected proteins")
legend("topright",row.names(protexp)[400:404], 
       lty = c(1,2,3,4,5), bty = "n", cex = .8)
dev.off()

## figure 11.7

tiff("figure11_7.tif",width=5,height=5,units="in",res=1200)
par(mai=c(0.4,0.8,0.4,0.1))
heatmap(protexp[350:404,],col=gray(100:1/100))
dev.off()