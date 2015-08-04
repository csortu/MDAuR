library(Cairo)

qpcr.data<-read.table("qPCR_intensities.txt")

CairoTIFF("figure5_2.tif",width=5,height=5,units="in",res=1200)
plot(qpcr.data,type="l",main="Change of intensity in a RT-PCR reaction")
dev.off()

qpcr.data2 <- read.table("dilutions.txt")

#CairoTIFF("/home/ortutay/Dropbox/Misc/BDAR_book/Wiley/Manuscript/05_Chapter_5/Ch_5_Figures/figure5_3.tif",width=5,height=5,units="in",res=1200)

tiff("figure5_3.tif",width=5,height=5,units="in",res=1200)
matplot(qpcr.data2,t="l",lty=c(rep(c(1,2,4),each=3),3,3),main="RT-PCR diluted samples", xlab="Cycles", ylab="Fluorescence",lwd=1,col="black")
legend("topleft", legend=c("Undiluted", "5 x diluted","25 x diluted","125 x diluted"),col=1,lty=c(1,2,4,3),lwd=1)
dev.off()

#Color version figure:
#matplot(qpcr.data2,t="l",col=c(rep(1:3,each=3),4,4),lty=1,main="RT-PCR diluted samples", xlab="Cycles", ylab="Fluorescence")
#legend("topleft", legend=c("Undiluted", "5 x diluted","25 x diluted","125 x diluted"),lty=1,col=1:4)

library(qpcR)
undil<-pcrfit(cbind(1:40,qpcr.data2),fluo=2:4)
tiff("figure5_4.tif",width=5,height=5,units="in",res=1200)
efficiency(undil)
dev.off()


dil5<-pcrfit(cbind(1:40,qpcr.data2),fluo=5:7)
dil25<-pcrfit(cbind(1:40,qpcr.data2),fluo=8:10)
dil125<-pcrfit(cbind(1:40,qpcr.data2),fluo=11:12)
undil.eff<-efficiency(undil,plot=F)
dil5.eff<-efficiency(dil5,plot=F)
dil25.eff<-efficiency(dil25,plot=F)
dil125.eff<-efficiency(dil125,plot=F)
threshold<-mean(c(undil.eff$fluo,dil5.eff$fluo,dil25.eff$fluo,dil125.eff$fluo))
Ct.undil <- efficiency(undil,plot=F,thres=threshold)$cpT
Ct.dil5 <- efficiency(dil5,plot=F,thres=threshold)$cpT
Ct.dil25 <- efficiency(dil25,plot=F,thres=threshold)$cpT
Ct.dil125 <- efficiency(dil125,plot=F,thres=threshold)$cpT




tiff("figure5_5.tif",width=5,height=5,units="in",res=1200)

matplot(qpcr.data2,t="l",col=1,lty=c(rep(1:3,each=3),4,4),main="RT-PCR diluted samples",xlab="Cycles",ylab="Fluorescence")
legend("topleft",legend=c("Undiluted","5 x diluted","25 x diluted","125 x diluted"),lty=1:4,col=1)
abline(h=threshold,lty=2)
lines(c(Ct.undil,Ct.undil),c(threshold,-300),lty=5,col=1)
lines(c(Ct.dil5,Ct.dil5),c(threshold,-300),lty=5,col=1)
lines(c(Ct.dil25,Ct.dil25),c(threshold,-300),lty=5,col=1)
lines(c(Ct.dil125,Ct.dil125),c(threshold,-300),lty=5,col=1)
text(10,1800,labels=paste("Threshold =",format(threshold,digits=0)))
text(c(Ct.undil,Ct.dil5,Ct.dil25,Ct.dil125),2500,labels=c(paste("Ct =",Ct.undil),paste("Ct =",Ct.dil5),paste("Ct =",Ct.dil25),paste("Ct =",Ct.dil125)),srt=90,col=1)
dev.off()

abs.dat <- read.table("abs_quant_data.txt")
abs.mod<-modlist(abs.dat,fluo=2:31)
abs.pars<-pcrbatch(abs.mod,group=rep(1:10,each=3),names="first")
threshold <- mean(as.numeric(abs.pars[abs.pars$Vars =="sig.fluo",-1]))
abs.pars.thres<-pcrbatch(abs.mod,group=rep(1:10,each=3),thres=threshold,names="first")
flu.data <- data.frame(Ct=as.numeric(abs.pars.thres[11,-1]), Load=c(10^(10:5), rep(NA,4)),row.names = c(paste("RefDil", rep(0:5),sep=""), paste("Sample",rep(1:4),sep="")))
flu.data$logLoad<-log10(flu.data$Load)

flu.mod<-lm(flu.data$logLoad ~ flu.data$Ct,subset=1:6)






tiff("figure5_6.tif",width=5,height=5,units="in",res=1200)

plot(flu.data$logLoad ~ flu.data$Ct,main="Viral load in influenza samples",xlab="Ct",ylab="log(Virus particles)")

text(flu.data$Ct[1:4],flu.data$logLoad[1:4]-0.7,labels=paste("RefDil",rep(0:3),sep=""),srt=90)
text(flu.data$Ct[5]-1.7,flu.data$logLoad[5]-0.5,labels="RefDil4",srt=45)
text(flu.data$Ct[6]-2.2,flu.data$logLoad[6],labels="RefDil5")

abline(flu.mod,col="black",lty=2)

points(flu.data[7:10,"Ct"],predict.lm(flu.mod,flu.data)[7:10],pch=5)
text(flu.data[7:9,"Ct"],predict.lm(flu.mod,flu.data)[7:9]+0.7,col="black",labels=paste("Sample",rep(1:3),sep=""),srt=90)
text(flu.data[10,"Ct"]+2.7,predict.lm(flu.mod,flu.data)[10],col="black",labels="Sample4")

dev.off()


library(RDML)
naive.rdml <- RDML$new("qPCR_data_naive_18S_IL2.rdml")
naive.info <- naive.rdml$AsTable()

sel.info <- naive.info[grepl("(koA 4 h)|(wt 4 h)",naive.info$sample),]
sel.info <- sel.info[order(sel.info$target,sel.info$sample),]
sel.dat <- naive.rdml$GetFData(sel.info,data.type="adp")


library(qpcR)

sel.mod <- modlist(sel.dat,fluo=2:13)
sel.pars <- pcrbatch(sel.mod)
replicates <- c("rs", "rs", "rs", "rc", "rc", "rc", "gs","gs", "gs", "gc", "gc", "gc")

tiff("figure5_7_col.tif",width=5,height=5,units="in",res=1200)

sel.ratio <- ratiocalc(sel.pars,group=replicates,type.eff = "mean.single")
mtext(c("A","B","C"),side=2,outer=F,at=c(1180,740,280),line=-1,las=1)

dev.off()

############################################################

library(RDML)
naive.rdml <- RDML$new("qPCR_data_naive_18S_IL2.rdml")
naive.info <- naive.rdml$AsTable()

all.info <- naive.info[grep("(24 h)|(48 h)",naive.info$sample,invert=T),]
all.info$sample <- sub(" 4 h"," 04 h",all.info$sample)
all.info <- all.info[order(all.info$target,all.info$sample),]
all.info[,5:7]
all.dat <- naive.rdml$GetFData(all.info,data.type="adp")

library(qpcR)

all.mod <- modlist(all.dat,fluo=2:55)
all.pars <- pcrbatch(all.mod)

groups2 <- paste(c(rep("r1",27), rep("g1",27)), rep(c(paste("s",1:6,sep=""), paste("c",1:3,sep="")), each=3), sep="")

all.ratio2 <- ratiobatch(all.pars,group=groups2,combs="same")

all.res2<-all.ratio2$resDat[,c(1,5,9,10,14,18)]
dimnames(all.res2)[[2]] <- c("KO A  0h","KO A  4h","KO A 12h","KO B  0h","KO B  4h","KO B 12h")

tiff("figure5_8.tif",width=5,height=5,units="in",res=1200)

plot(1:3-0.1,all.res2[7,1:3],ylim=c(min(all.res2[11,]),max(all.res2[12,])),main="IL-2 activation in KO mice",xaxt="n",xlab="Time after activation",ylab="IL-2 expression in KO/WT",xlim=c(0.5,3.5))
axis(1,at=1:3,labels=c("0h","4h","12h"),tick=F)
arrows(1:3-0.1,all.res2[11,1:3],1:3-0.1,all.res2[12,1:3], length=0.05, angle=90, code=3)

points(1:3+0.1,all.res2[7,4:6],pch=18)
arrows(1:3+0.1,all.res2[11,4:6],1:3+0.1,all.res2[12,4:6], length=0.05, angle=90, code=3,lty=2)

legend(2.6,1.7,c("KO A","KO B"),pch=c(1,18),lty=c(1,2))

abline(h=1,lty=3)

dev.off()

############################################################

library(RDML)

naive.rdml <- RDML$new("qPCR_data_naive_18S_IL2.rdml")
naive.info <- naive.rdml$AsTable()
naive <- naive.rdml$GetFData(naive.info,data.type="mdp")

library(qpcR)

tiff("figure5_9_color.tif",width=5,height=5,units="in",res=1200)

melt <- meltcurve(naive[,c(1,2)],cut.Area =10)

dev.off()

