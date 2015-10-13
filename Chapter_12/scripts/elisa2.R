library(readxl)

xls.file <- "090714_data.xlsx"
excel_sheets(xls.file)

ch450 <- read_excel(xls.file,sheet=1)
ch570 <- read_excel(xls.file,sheet=2)
samples <- read_excel(xls.file,sheet=3)

ch450[,1] <- NULL
ch570[,1] <- NULL
samples[,1] <- NULL

elisa.data <- data.frame(ch450=unlist(ch450),ch570=unlist(ch570),sample=unlist(samples))

summary(elisa.data$ch450)
summary(elisa.data$ch570)

image(1:12,1:8,t(matrix(elisa.data$ch570,nrow=8)),xlab="",ylab="",ylim=c(8.5,0.5),main="Channel 570")
text(rep(1:12,each=8),rep(1:8,12),labels=matrix(elisa.data$sample,nrow=8),pos=3)
text(rep(1:12,each=8),rep(1:8,12),labels=matrix(elisa.data$ch570,nrow=8),pos=1)


image(1:12,1:8,t(matrix(elisa.data$ch450,nrow=8)),xlab="",ylab="",ylim=c(8.5,0.5),main="Channel 450")
text(rep(1:12,each=8),rep(1:8,12),labels=matrix(elisa.data$sample,nrow=8),pos=3)
text(rep(1:12,each=8),rep(1:8,12),labels=matrix(elisa.data$ch450,nrow=8),pos=1)

#elisa.data$platecol <- rep(1:12,each=8)
#elisa.data$platerow <- as.factor(rep(LETTERS[1:8],12))

elisa.data$od <- elisa.data$ch450 - elisa.data$ch570

image(1:12,1:8,t(matrix(elisa.data$od,nrow=8)),xlab="",ylab="",ylim=c(8.5,0.5),main="OD",col=colorRampPalette(c("white","yellow"))(20))
text(rep(1:12,each=8),rep(1:8,12),labels=matrix(elisa.data$sample,nrow=8),pos=3)
text(rep(1:12,each=8),rep(1:8,12),labels=matrix(elisa.data$od,nrow=8),pos=1)

## fit linear model to standards

elisa.data$sample[grep("Std",elisa.data$sample)]


elisa.standard <- data.frame(lab=elisa.data$sample[grepl("Std",elisa.data$sample)],od=elisa.data$od[grepl("Std",elisa.data$sample)])

#elisa.standard$conc<-sub(" pg/ml","",sub("Std ","",elisa.standard$lab))
#conc.pat <- "(\\d)+"

library(stringr)

conc.pat <- "(\\d+\\.?\\d+)"
  
elisa.standard$conc <- as.numeric(str_extract(elisa.standard$lab,conc.pat))

elisa.st.lm <- lm(conc ~ od - 1, data=elisa.standard)

plot(elisa.standard$conc ~ elisa.standard$od,ylab="IL2 concentration",xlab="OD",main="ELISA standard",sub="linear model")
abline(elisa.st.lm)

plot(elisa.st.lm)

## 4pl model

library(drc)

elisa.st.4pl <- drm(formula = conc ~ od, data = elisa.standard, fct = LL.4())

plot(elisa.standard$conc ~ elisa.standard$od,ylab="IL2 concentration",xlab="OD",main="ELISA standard",sub="Four-parameter logistic model")
plot(elisa.st.4pl,add=T)
abline(elisa.st.lm,lty=2)

elisa.data$conc.lm <- predict(elisa.st.lm,data.frame(od=elisa.data$od))
elisa.data$conc.4pl <- predict(elisa.st.4pl,data.frame(od=elisa.data$od))

## evaluation

sum(elisa.data$conc.lm<0)
sum(elisa.data$conc.4pl<0)


## too many extrapolation! repeat the experiment with 

elisa.data$color <- "green"
elisa.data$color[elisa.data$conc.4pl>200] <- "red"

image(1:12,1:8,t(matrix(as.numeric(as.factor(elisa.data$color)),nrow=8)),xlab="",ylab="",ylim=c(8.5,0.5),main="OD",col=c("green","red"))
text(rep(1:12,each=8),rep(1:8,12),labels=matrix(elisa.data$sample,nrow=8),pos=3)
text(rep(1:12,each=8),rep(1:8,12),labels=format(matrix(elisa.data$conc.4pl,nrow=8),digits=2),pos=1)
#text(rep(1:12,each=8),rep(1:8,12),labels=format(matrix(elisa.data$conc.4pl,nrow=8),digits=2),pos=4)
#text(rep(1:12,each=8),rep(1:8,12),labels=format(matrix(elisa.data$conc.lm,nrow=8),digits=0),pos=2)


