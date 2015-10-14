library(readxl)
xls.file <- "090714_data.xlsx"
ch450 <- read_excel(xls.file,sheet=1)
ch570 <- read_excel(xls.file,sheet=2)
samples <- read_excel(xls.file,sheet=3)
ch450[,1] <- NULL
ch570[,1] <- NULL
samples[,1] <- NULL

elisa.data <- data.frame(ch450=unlist(ch450), 
                         ch570=unlist(ch570), 
                         sample=unlist(samples))

## figure 12.2

tiff("figure12_2.tif",width=5,height=5,units="in",res=1200)
image(1:12,1:8,
      t(matrix(elisa.data$ch570,nrow=8)),
      xlab="",
      ylab="",ylim=c(8.5,0.5),
      main="Channel 570",
      col = gray.colors(12))
text(rep(1:12,each=8),
     rep(1:8,12),
     labels=matrix(elisa.data$sample,nrow=8),
     pos=3,cex=0.35)
text(rep(1:12,each=8),
     rep(1:8,12),
     labels=matrix(elisa.data$ch570,nrow=8),
     pos=1,cex=0.4)
dev.off()

elisa.data$od <- elisa.data$ch450 - elisa.data$ch570
elisa.data$sample[grep("Std",elisa.data$sample)]
elisa.standard <- data.frame(lab=elisa.data$sample[grepl("Std",elisa.data$sample)], od=elisa.data$od[grepl("Std",elisa.data$sample)])
elisa.standard$conc <- c(200, 100, 50, 25, 12.5, 6.25,200, 100, 50, 25, 12.5, 6.25)

elisa.st.lm <- lm(conc ~ od - 1, data=elisa.standard)

## figure 12.3
tiff("figure12_3.tif",width=5,height=5,units="in",res=1200)
par(mfrow=c(2,1),mai=c(0.77,0.82,0.52,0.12))
plot(elisa.standard$conc ~ elisa.standard$od,ylab="IL2 concentration",xlab="OD",main="ELISA standard",sub="linear model")
abline(elisa.st.lm)
mtext("A", side = 3, line = 0, adj = 0, cex = 1.3)
plot(elisa.st.lm,which=1,caption=NA)
mtext("B", side = 3, line = 0, adj = 0, cex = 1.3)
dev.off()

library(drc)
elisa.st.4pl <- drm(formula = conc ~ od, data = elisa.standard, fct = LL.4())

## figure 12.4
tiff("figure12_4.tif",width=5,height=5,units="in",res=1200)
par(mfrow=c(1,1),mai=c(0.77,0.82,0.52,0.12))
plot(elisa.standard$conc ~ elisa.standard$od,
     ylab="IL2 concentration",xlab="OD",
     main="ELISA standard")
plot(elisa.st.4pl,add=T)
abline(elisa.st.lm,lty=2)
legend("topleft",
       c("linear model","Four-parameter logistic model"),
       bty="n",
       lty=c(1,2))
dev.off()

## figure 12.5

library(readxl)
xls.file <- "240214_data.xlsx"
ch450 <- read_excel(xls.file,sheet=1)
ch570 <- read_excel(xls.file,sheet=2)
samples <- read_excel(xls.file,sheet=3)
ch450[,1] <- NULL
ch570[,1] <- NULL
samples[,1] <- NULL
elisa.data <- data.frame(ch450=unlist(ch450),ch570=unlist(ch570),sample=unlist(samples))
elisa.data$od <- elisa.data$ch450 - elisa.data$ch570
elisa.standard <- data.frame(lab=elisa.data$sample[grepl("IL-2",elisa.data$sample)],od=elisa.data$od[grepl("IL-2",elisa.data$sample)])
library(stringr)
conc.pat <- "(\\d+\\.?\\d+)"
elisa.standard$conc <- as.numeric(str_extract(elisa.standard$lab,conc.pat))
library(drc)
elisa.st.4pl <- drm(formula = conc ~ od, data = elisa.standard, fct = LL.4())
elisa.data$conc.4pl <- predict(elisa.st.4pl,data.frame(od=elisa.data$od))

elisa.data$animal[grep("IL-2",elisa.data$sample,invert=T)] <- unlist(strsplit(as.character(elisa.data$sample[grep("IL-2", elisa.data$sample,invert=T)]), " "))[seq(from=1, by=4, length.out=80)]
elisa.data$dilution[grep("IL-2",elisa.data$sample,invert=T)] <- unlist(strsplit(as.character(elisa.data$sample[grep("IL-2", elisa.data$sample,invert=T)]), " "))[seq(from=2, by=4, length.out=80)]
elisa.data$treatment[grep("IL-2",elisa.data$sample,invert=T)] <- unlist(strsplit(as.character(elisa.data$sample[grep("IL-2", elisa.data$sample,invert=T)]), " "))[seq(from=3, by=4, length.out=80)]
elisa.data$time[grep("IL-2",elisa.data$sample,invert=T)] <- unlist(strsplit(as.character(elisa.data$sample[grep("IL-2", elisa.data$sample,invert=T)]), " "))[seq(from=4, by=4, length.out=80)]
elisa.data$genotype <- substr(elisa.data$animal,1,2)

measured.data <- elisa.data[elisa.data$time=="4h"&!is.na(elisa.data$animal),c(10:6,5)]
measured.data <- rbind(measured.data,elisa.data[elisa.data$time=="24h"&elisa.data$dilution=="1/100"&!is.na(elisa.data$animal),c(10:6,5)])
measured.data <- rbind(measured.data,elisa.data[elisa.data$time=="48h"&elisa.data$dilution=="1/1000"&!is.na(elisa.data$animal),c(10:6,5)])

measured.data$undil.conc <- measured.data$conc.4pl
measured.data$undil.conc[measured.data$dilution=="1/100"] <- measured.data$undil.conc[measured.data$dilution=="1/100"] * 100
measured.data$undil.conc[measured.data$dilution=="1/1000"] <- measured.data$undil.conc[measured.data$dilution=="1/1000"] * 1000

measured.data$time <- as.factor(measured.data$time)
measured.data$time <- relevel(measured.data$time,"4h")

tiff("figure12_5.tif",width=5,height=5,units="in",res=1200)
boxplot(undil.conc ~ time,data=measured.data,log="y",
        main="IL-2 production of T-cells\ntime series",
        ylab="IL-2 [pg/ml]")
dev.off()

## figure 12.6

tiff("figure12_6.tif",width=5,height=5,units="in",res=1200)
par(xaxt="n")
boxplot(undil.conc ~ time + genotype, 
        data=measured.data[measured.data$treatment=="Tr1",], 
        log="y", 
        main="IL-2 production of T-cells in FURIN WT/KO mice", 
        ylab="IL-2 [pg/ml]")
lablist <- rep(c("4h","24h","48h"),2)
text(1:6, 4, labels = lablist, srt = 0, pos = 1, xpd = TRUE)
lablist <- rep(c("KO","WT"),each=3)
text(1:6, 8, labels = lablist, srt = 0, pos = 1, xpd = TRUE)
dev.off() 
