library(readxl)

xls.file <- "240214_data.xlsx"
excel_sheets(xls.file)

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

## elisa.data$color <- "green"
## elisa.data$color[elisa.data$conc.4pl>400] <- "red"

## image(1:12,1:8,t(matrix(as.numeric(as.factor(elisa.data$color)),nrow=8)),xlab="",ylab="",ylim=c(8.5,0.5),main="OD",col=c("green","red"))
## text(rep(1:12,each=8),rep(1:8,12),labels=matrix(elisa.data$sample,nrow=8),pos=3)
## text(rep(1:12,each=8),rep(1:8,12),labels=format(matrix(elisa.data$conc.4pl,nrow=8),digits=2),pos=1)

## comparative analysis


elisa.data[grep("IL-2",elisa.data$sample,invert=T),]

elisa.data$animal[grep("IL-2",elisa.data$sample,invert=T)] <- unlist(strsplit(as.character(elisa.data$sample[grep("IL-2",elisa.data$sample,invert=T)])," "))[seq(from=1,by=4,length.out=80)]

elisa.data$dilution[grep("IL-2",elisa.data$sample,invert=T)] <- unlist(strsplit(as.character(elisa.data$sample[grep("IL-2",elisa.data$sample,invert=T)])," "))[seq(from=2,by=4,length.out=80)]

elisa.data$treatment[grep("IL-2",elisa.data$sample,invert=T)] <- unlist(strsplit(as.character(elisa.data$sample[grep("IL-2",elisa.data$sample,invert=T)])," "))[seq(from=3,by=4,length.out=80)]

elisa.data$time[grep("IL-2",elisa.data$sample,invert=T)] <- unlist(strsplit(as.character(elisa.data$sample[grep("IL-2",elisa.data$sample,invert=T)])," "))[seq(from=4,by=4,length.out=80)]

elisa.data$genotype <- substr(elisa.data$animal,1,2)

## check proper dilutions

elisa.data[elisa.data$conc.4pl<400&!is.na(elisa.data$animal),c(10:7,5)]

##  4h -> 1
## 24h -> 1/100
## 48h -> 1/1000

elisa.data[elisa.data$time=="24h"&elisa.data$dilution=="1/100"&!is.na(elisa.data$animal),c(10:6,5)]

measured.data <- elisa.data[elisa.data$time=="4h"&!is.na(elisa.data$animal),c(10:6,5)]

measured.data <- rbind(measured.data,elisa.data[elisa.data$time=="24h"&elisa.data$dilution=="1/100"&!is.na(elisa.data$animal),c(10:6,5)])

measured.data <- rbind(measured.data,elisa.data[elisa.data$time=="48h"&elisa.data$dilution=="1/1000"&!is.na(elisa.data$animal),c(10:6,5)])

measured.data$undil.conc <- measured.data$conc.4pl

measured.data$undil.conc[measured.data$dilution=="1/100"] <- measured.data$undil.conc[measured.data$dilution=="1/100"] * 100

measured.data$undil.conc[measured.data$dilution=="1/1000"] <- measured.data$undil.conc[measured.data$dilution=="1/1000"] * 1000

str(measured.data)
measured.data$time <- as.factor(measured.data$time)

measured.data <- as.data.frame(unclass(measured.data)) # excellent shortcut
str(measured.data)

levels(measured.data$time)

measured.data$time <- relevel(measured.data$time,"4h")

boxplot(undil.conc ~ time,data=measured.data)
boxplot(undil.conc ~ time,data=measured.data,log="y")

boxplot(undil.conc ~ time + genotype,data=measured.data,log="y")
boxplot(undil.conc ~ time + treatment,data=measured.data,log="y")

boxplot(undil.conc ~ time + genotype + treatment,data=measured.data,log="y")
boxplot(undil.conc ~ time + treatment + genotype,data=measured.data,log="y")

boxplot(undil.conc ~ time + genotype,data=measured.data[measured.data$treatment=="Tr1",],log="y",main="IL-2 production of T-cells in FURIN KO mice",ylab="IL-2 [pg/ml]")
boxplot(undil.conc ~ time + genotype,data=measured.data[measured.data$treatment=="10+10",],log="y",main="IL-2 production of T-cells in FURIN KO mice",ylab="IL-2 [pg/ml]")
