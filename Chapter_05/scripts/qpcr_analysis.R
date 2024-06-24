qpcr.data<-read.table("qPCR_intensities.txt")
plot(qpcr.data,type="l")

qpcr.data2 <- read.table("dilutions.txt")
matplot(qpcr.data2,t="l",col=c(rep(1:3,each=3),4,4),lty=1,main="RT-PCR diluted samples",xlab="Cycles",ylab="Fluorescence")
legend("topleft",legend=c("Undiluted","5 x diluted","25 x diluted","125 x diluted"),lty=1,col=1:4)

# dCt demo
matplot(qpcr.data2[,c(1,7,11)],t="l",col=2:4,lty=1)
abline(v=27,lty=2)
qpcr.data2[,7]/qpcr.data2[,11]

matplot(qpcr.data2[,c(1,7,11)],t="l",col=2:4,lty=1)
abline(h=2000,lty=2)
lines(c(24,24),c(2333,-300),lty=2,col="red")
lines(c(28,28),c(2130,-300),lty=2,col="green")
lines(c(31,31),c(2110,-300),lty=2,col="blue")
text(c(24,28,31),3000,labels=c("Ct=24","Ct=28","Ct=31"),srt=90,col=2:4)


library(qpcR)

undil<-pcrfit(cbind(1:40,qpcr.data2),fluo=2:4)
dil5<-pcrfit(cbind(1:40,qpcr.data2),fluo=5:7)
dil25<-pcrfit(cbind(1:40,qpcr.data2),fluo=8:10)
dil125<-pcrfit(cbind(1:40,qpcr.data2),fluo=11:12)

plot(undil)

efficiency(undil)

undil.eff<-efficiency(undil,plot=F)
dil5.eff<-efficiency(dil5,plot=F)
dil25.eff<-efficiency(dil25,plot=F)
dil125.eff<-efficiency(dil125,plot=F)

threshold<-mean(c(undil.eff$fluo,dil5.eff$fluo,dil25.eff$fluo,dil125.eff$fluo))

Ct.undil <- efficiency(undil,plot=F,thres=threshold)$cpT
Ct.dil5 <- efficiency(dil5,plot=F,thres=threshold)$cpT
Ct.dil25 <- efficiency(dil25,plot=F,thres=threshold)$cpT
Ct.dil125 <- efficiency(dil125,plot=F,thres=threshold)$cpT

matplot(qpcr.data2,t="l",col=c(rep(1:3,each=3),4,4),lty=1,main="RT-PCR diluted samples",xlab="Cycles",ylab="Fluorescence")
legend("topleft",legend=c("Undiluted","5 x diluted","25 x diluted","125 x diluted"),lty=1,col=1:4)
abline(h=threshold,lty=2)
lines(c(Ct.undil,Ct.undil),c(threshold,-300),lty=2,col="black")
lines(c(Ct.dil5,Ct.dil5),c(threshold,-300),lty=2,col="red")
lines(c(Ct.dil25,Ct.dil25),c(threshold,-300),lty=2,col="green")
lines(c(Ct.dil125,Ct.dil125),c(threshold,-300),lty=2,col="blue")
text(10,1800,labels=paste("Threshold =",threshold))
text(c(Ct.undil,Ct.dil5,Ct.dil25,Ct.dil125),2500,labels=c(paste("Ct =",Ct.undil),paste("Ct =",Ct.dil5),paste("Ct =",Ct.dil25),paste("Ct =",Ct.dil125)),srt=90,col=1:4)

((undil.eff$eff+dil5.eff$eff)/2)^(Ct.dil5-Ct.undil)







library(RDML)

## This is for RDML version before ver 0.8-4 (2015-04-19)
## calib <- RDML("calibration.rdml")
## plot(1:40,calib$qPCR$"18s RNA"$unkn[,14],type="l",main="KO")
## lines(1:40,calib$qPCR$"18s RNA"$unkn[,17],col="red")
## lines(1:40,calib$qPCR$"18s RNA"$unkn[,20],col="blue")
## lines(1:40,calib$qPCR$"18s RNA"$unkn[,23],col="green")

## naive <- RDML("qPCR_data_naive_18S_IL2.rdml")
## matplot(naive$qPCR$"18s"$unkn[,2:40],type="l")
## names(naive)



naive.rdml <- RDML$new("qPCR_data_naive_18S_IL2.rdml")
naive.info <- naive.rdml$AsTable()
str(naive.info)
naive.info$fdata.name
naive.info$sample
naive.info$target
naive.info$position

# # This is for version before ver 0.9-3
# # parameter data.type was renamed to dp.type
# naive.rdml$GetFData(naive.info,data.type="mdp")
# matplot(naive.rdml$GetFData(naive.info[naive.info$target=="18s",],data.type="adp"),type="l")

naive.rdml$GetFData(naive.info,dp.type="mdp")
matplot(naive.rdml$GetFData(naive.info[naive.info$target=="18s",],dp.type="adp"),type="l")



## model18S<-modlist(calib$qPCR$"18s RNA"$unkn,fluo=2:24)
## modelFgf2<-modlist(calib$qPCR$"Fgf-2"$unkn,fluo=2:24)
## modelIl2<-modlist(calib$qPCR$"IL-2"$unkn,fluo=2:24)
## plot(model18S,which="3D",col=rep(1:8, each = 3))
## rgl.close()
## plot(modelFgf2,which="3D",col=rep(1:8, each = 3))
## rgl.close()
## plot(modelIl2,which="3D",col=rep(1:8, each = 3))
## rgl.close()
## naivemodel18s <- modlist(naive$qPCR$"18s"$unkn,fluo=2:40)
## naivemodelIl2 <- modlist(naive$qPCR$"IL-2"$unkn,fluo=2:40)
## plot(naivemodelIl2,which="3D")
## rgl.close()


#pcrbatch
#ratiocalc


# Absolute quantification

abs.dat <- read.table("abs_quant_data.txt")

library(qpcR)

abs.mod<-modlist(abs.dat,fluo=2:31)
plot(abs.mod,col=rep(1:10,each=3))


class(abs.mod) <- "modlist"

abs.pars<-pcrbatch(abs.mod,
                   group=rep(1:10,each=3),
                   names="first")


abs.pars[c(8,15),]
threshold <- mean(as.numeric(abs.pars[abs.pars$Vars =="sig.fluo",-1]))


abs.pars.thres<-pcrbatch(abs.mod,
                         group=rep(1:10,each=3),
                         thres=threshold,
                         names="first")

abs.pars.thres[c(8,11),]

# Fit linear model and find the sample concentrations

flu.data <- data.frame(Ct=as.numeric(abs.pars.thres[11,-1]),Load=c(10^(10:5),rep(NA,4)),row.names = c(paste("RefDil",rep(0:5),sep=""),paste("Sample",rep(1:4),sep="")))

flu.data$logLoad<-log10(flu.data$Load)

par(mfrow=c(2,1))
plot(flu.data$Load ~ flu.data$Ct)
plot(flu.data$logLoad ~ flu.data$Ct)
par(mfrow=c(1,1))
  
flu.mod<-lm(flu.data$logLoad ~ flu.data$Ct,subset=1:6)

plot(flu.data$logLoad ~ flu.data$Ct)
abline(flu.mod,col="red")

predict.lm(flu.mod,flu.data)
predict.lm(flu.mod,flu.data)[7:10]

flu.data$logLoad[7:10]<-predict.lm(flu.mod,flu.data)[7:10]
flu.data$Load[7:10] <- 10^flu.data$logLoad[7:10]



plot(flu.data$logLoad ~ flu.data$Ct,main="Viral load in influenza samples",xlab="Ct",ylab="log(Virus particles)")
text(flu.data$Ct[1:4],flu.data$logLoad[1:4]-0.5,labels=paste("RefDil",rep(0:3),sep=""),srt=90)
text(flu.data$Ct[5:6]-2,flu.data$logLoad[5:6],labels=paste("RefDil",4:5,sep=""))
abline(flu.mod,col="red")
points(flu.data[7:10,"Ct"],predict.lm(flu.mod,flu.data)[7:10],col="blue")
text(flu.data[7:10,"Ct"],predict.lm(flu.mod,flu.data)[7:10]+0.5,col="blue",labels=paste("Sample",rep(1:4),sep=""),srt=90)

for (i in 7:10){
  lines(c(0,flu.data$Ct[i],flu.data$Ct[i]),c(flu.data$logLoad[i],flu.data$logLoad[i],0),lty=2,col="blue")
}

# Quick absolute quantification

abs.ref.mod <- modlist(abs.dat,fluo=2:19)
abs.sampl.mod <- modlist(abs.dat,fluo=20:31)

par(mfrow=c(2,1))
plot(abs.ref.mod,col=rep(1:6,each=3))
plot(abs.sampl.mod,col=rep(1:4,each=3))
par(mfrow=c(1,1))

preds <- calib(abs.ref.mod,predcurve = abs.sampl.mod,dil=10^(10:5), group = gl(6,3))

for(i in 1:4){
  print(10^mean(preds$predconc[,(i*3-2):(i*3)]))
}

boxplot(preds$predconc,main="Predictions of sample concentrations",xlab="Samples",ylab="log Viral load")


# Relative quantification with ratiocalc

library(RDML)

##  # old RDML version
## naive <- RDML("qPCR_data_naive_18S_IL2.rdml")
## # two conditions only
## sel.dat<-cbind(naive$qPCR$"18s"$unkn[,c(1,7,19,38,8,20,39)],naive$qPCR$"IL-2"$unkn[,c(2,16,28,3,17,29)])
## names(sel.dat)<-c("Cycles","18s.wt.1","18s.wt.2","18s.wt.3","18s.ko.1","18s.ko.2","18s.ko.3","Il2.wt.1","Il2.wt.2","Il2.wt.3","Il2.ko.1","Il2.ko.2","Il2.ko.3")


naive.rdml <- RDML$new("qPCR_data_naive_18S_IL2.rdml")
naive.info <- naive.rdml$AsTable()

#naive.rdml$GetFData(naive.info[grepl("(koA)|(wt)",naive.info$sample),],data.type="adp")
#naive.info[grepl("(koA 4 h)|(wt 4 h)",naive.info$sample),]

sel.info <- naive.info[grep("(koA 4 h)|(wt 4 h)",naive.info$sample),]
sel.info <- sel.info[order(sel.info$target,sel.info$sample),]

# # This is for version before ver 0.9-3
# # parameter data.type was renamed to dp.type
# sel.dat <- naive.rdml$GetFData(sel.info,data.type="adp")
sel.dat <- naive.rdml$GetFData(sel.info,dp.type="adp")
colnames(sel.dat)


library(qpcR)

# sel.mod <- modlist(sel.dat,fluo=2:13)

#sel.pars <- pcrbatch(sel.mod,group=rep(1:4,each=3),names="first")
#sel.ratio <- ratiocalc(sel.pars,group=c("rc","rs","gc","gs"))

#sel.pars <- pcrbatch(sel.mod)


# Latest version of pcrbatch() has issues with handling data 
# with multiple class labels. 
# Here we are falling back to use the data frame

class(sel.dat) <- "data.frame"
names(sel.dat)[1] <- 'Cycles'
sel.pars <- pcrbatch(sel.dat,
                     cyc = 1,
                     fluo=2:13)

## replicates <- c("rc","rc","rc","rs","rs","rs","gc","gc","gc","gs","gs","gs")

replicates <- c("rs", "rs", "rs", "rc", "rc", "rc", "gs","gs", "gs", "gc", "gc", "gc")
sel.ratio <- ratiocalc(sel.pars,group=replicates,type.eff = "mean.single")

summary(sel.ratio)
sel.ratio$summary

# all data with two genes, two conditions: wt and koA + koB (replicates) at different time points: 0h, 4h, 12h three technical replicates for most conditions WT 0h is used as universal calibrator condition

all.info <- naive.info[grep("(24 h)|(48 h)",naive.info$sample,invert=T),]
all.info$sample <- sub(" 4 h"," 04 h",all.info$sample)
all.info <- all.info[order(all.info$target,all.info$sample),]
all.info[,5:7]

# # This is for version before ver 0.9-3
# # parameter data.type was renamed to dp.type
# all.dat <- naive.rdml$GetFData(all.info,data.type="adp")
all.dat <- naive.rdml$GetFData(all.info,dp.type="adp")
colnames(all.dat)


## old RDML version
## all.dat <- cbind(naive$qPCR$"18s"$unkn[,c(1,         # Cycle count
##                                           5,17,29,   # WT    0h
##                                           7,19,38,   # WT    4h
##                                           10,22,31,  # WT   12h
##                                           6,18,30,   # KO A  0h
##                                           8,20,39,   # KO A  4h
##                                           11,23,32,  # KO A 12h
##                                           2,3,4,     # KO B  0h
##                                           9,21,14,   # KO B  4h
##                                           12,24,33   # KO B 12h
##                                           )],
##                  naive$qPCR$"IL-2"$unkn[,c(5,14,26,  # WT    0h
##                                            2,16,28,  # WT    4h
##                                            7,19,31,  # WT   12h
##                                            6,15,27,  # KO A  0h
##                                            3,17,29,  # KO A  4h
##                                            8,20,32,  # KO A 12h
##                                            38,39,40, # KO B  0h
##                                            4,18,30,  # KO B  4h
##                                            9,21,33   # KO B 12h
##                                            )])

## names(all.dat) <- c("Cycle",paste(rep(c("18s","Il2"),each=27),".",rep(c("WT","KoA","KoB"),each=9),".",rep(c("0h","4h","12h"),each=3),".",1:3,sep=""))

# all.mod <- modlist(all.dat,fluo=2:55)
# all.pars <- pcrbatch(all.mod)

# Latest version of pcrbatch() has issues with handling data 
# with multiple class labels. 
# Here we are falling back to use the data frame
class(all.dat) <- "data.frame"
names(all.dat)[1] <- 'Cycles'

all.pars <- pcrbatch(all.dat,fluo=2:55)

##groups2 <- paste(c(rep("r1",27),rep("g1",27)),rep(c(paste("s",1:6,sep=""),paste("c",1:3,sep="")),each=3),sep="")

groups2 <- paste(c(rep("r1",27), rep("g1",27)), rep(c(paste("s",1:6,sep=""), paste("c",1:3,sep="")), each=3), sep="")


all.ratio2 <- ratiobatch(all.pars,group=groups2,combs="same")

all.res2<-all.ratio2$resDat[,c(1,5,9,10,14,18)]
dimnames(all.res2)[[2]] <- c("KO A  0h","KO A  4h","KO A 12h","KO B  0h","KO B  4h","KO B 12h")

dev.off()
plot(1:3-0.1,all.res2[7,1:3],ylim=c(min(all.res2[11,]),max(all.res2[12,])),main="IL-2 activation in KO mice",xaxt="n",xlab="Time after activation",ylab="IL-2 expression in KO/WT",xlim=c(0.5,3.5))
axis(1,at=1:3,labels=c("0h","4h","12h"),tick=F)
arrows(1:3-0.1,all.res2[11,1:3],1:3-0.1,all.res2[12,1:3], length=0.05, angle=90, code=3)

points(1:3+0.1,all.res2[7,4:6],col="blue")
arrows(1:3+0.1,all.res2[11,4:6],1:3+0.1,all.res2[12,4:6], length=0.05, angle=90, code=3,col="blue")

#legend(2.8,1.5,c("KO A","KO B"),text.col=c("black","blue"),lty=1,col=c("black","blue"))
legend(2.3,50,c("KO A","KO B"),text.col=c("black","blue"),lty=1,col=c("black","blue"))


# Meltcurve analysis

library(RDML)

## Old RDML version
## naive <- RDML("qPCR_data_naive_18S_IL2.rdml")
## melt <- meltcurve(naive$Melt$"18s"$unkn[,c(1,2)],cut.Area =10)

naive.rdml <- RDML$new("qPCR_data_naive_18S_IL2.rdml")
naive.info <- naive.rdml$AsTable()
# # This is for version before ver 0.9-3
# # parameter data.type was renamed to dp.type
# naive <- naive.rdml$GetFData(naive.info,data.type="mdp")
naive <- naive.rdml$GetFData(naive.info,dp.type="mdp")

library(qpcR)

melt <- meltcurve(naive[,c(1,2)],cut.Area =10)

head(melt[[1]])
melt[[1]]$Tm[1]


meltlist <- meltcurve(naive[,c(as.vector(rbind(rep(1,39),2:40)))],
                      plot = FALSE)

for (i in 1:length(meltlist)){
  print(meltlist[[i]]$Tm[ which(meltlist[[i]]$Area==max(meltlist[[i]]$Area,na.rm=T))])
}
