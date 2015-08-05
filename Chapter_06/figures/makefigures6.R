#datapath <- '/home/data/work/Courses/Microarray/GSE29797/'
setwd(datapath)
library(affy)
dat<-ReadAffy()
exp.des <- pData(dat)
exp.des$Genotype <- factor(rep(c("WT","dTsc1"),each=6))
exp.des$Stimulation <- factor(rep(c("0h","4h"),6))
exp.des
pData(dat) <- exp.des

library(Cairo)

CairoTIFF("figure6_2.tif",width=5,height=5,units="in",res=1200)
par(mfrow=c(2,2),mai=c(0.42,0.67,0.42,0.12))
image(dat[,1],main="Array image")
mtext("A", side = 3, line = -0.7, adj = -0.3, cex = 1.3)
hist(dat[,1],col=1,main="Raw intensities\nsingle sample")
mtext("B", side = 3, line = -0.7, adj = -0.3, cex = 1.3)
boxplot(exprs(dat),names=1:12,pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.3),xlab="samples",main="Raw intensities")
mtext("C", side = 3, line = -0.7, adj = -0.3, cex = 1.3)
boxplot(log2(exprs(dat)),names=1:12,pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.3),xlab="samples",main="Log intensities")
mtext("D", side = 3, line = -0.7, adj = -0.3, cex = 1.3)
dev.off()


#par(mfrow=c(1,1))

datrma<-rma(dat)

matexp<-exprs(datrma)

CairoTIFF("figure6_3.tif",width=5,height=5,units="in",res=1200)
par(mfrow=c(1,2),mai=c(0.52,0.77,0.52,0.32))
boxplot(log2(exprs(dat)),names=1:12,pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.3),main="Before normalization",ylab="log2(Intensity)",cex=0.7)
boxplot(matexp,names=1:12,pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.3),main="RMA normalized data",cex=0.7)
dev.off()


library(affycoretools)

CairoTIFF("figure6_4.tif",width=5,height=5,units="in",res=1200)

par(mfrow=c(1,2))
plotPCA(matexp, groups = as.numeric(pData(dat)[,2]), groupnames = levels(pData(dat)[,2]),col=c(1,1),pch=c(1,18),main="Principal Components\nWT vs dTsc1")
plotPCA(matexp, groups = as.numeric(pData(dat)[,3]), groupnames = levels(pData(dat)[,3]),col=c(1,1),pch=c(1,18),main="Principal Components\n0h vs 4h")

dev.off()

