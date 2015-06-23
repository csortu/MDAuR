library(Cairo)

qpcr.data<-read.table("qPCR_intensities.txt")

CairoTIFF("figure5_2.tif",width=5,height=5,units="in",res=1200)
plot(qpcr.data,type="l",main="Change of intensity in a RT-PCR reaction")
dev.off()

qpcr.data2 <- read.table("dilutions.txt")

#CairoTIFF("/home/ortutay/Dropbox/Misc/BDAR_book/Wiley/Manuscript/05_Chapter_5/Ch_5_Figures/figure5_3.tif",width=5,height=5,units="in",res=1200)

  tiff("/home/ortutay/Dropbox/Misc/BDAR_book/Wiley/Manuscript/05_Chapter_5/Ch_5_Figures/figure5_3.tif",width=5,height=5,units="in",res=1200)
  matplot(qpcr.data2,t="l",lty=c(rep(c(1,2,4),each=3),3,3),main="RT-PCR diluted samples", xlab="Cycles", ylab="Fluorescence",lwd=1,col="black")
  legend("topleft", legend=c("Undiluted", "5 x diluted","25 x diluted","125 x diluted"),col=1,lty=c(1,2,4,3),lwd=1)
  dev.off()

#Color version figure:
#matplot(qpcr.data2,t="l",col=c(rep(1:3,each=3),4,4),lty=1,main="RT-PCR diluted samples", xlab="Cycles", ylab="Fluorescence")
#legend("topleft", legend=c("Undiluted", "5 x diluted","25 x diluted","125 x diluted"),lty=1,col=1:4)

library(qpcR)
undil<-pcrfit(cbind(1:40,qpcr.data2),fluo=2:4)
efficiency(undil)
