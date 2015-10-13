library(drc)

stand<-data.frame(concentration=c(200,100,50,25,12.5,6.25),OD=c(2.1215,1.301,0.7,0.378,0.1905,0.104))

FourP<-drm(formula = concentration ~ OD, data = stand, fct = LL.4())
predict(FourP,data.frame(x=seq(0,1, by=0.1)  ))

#measured.data<-data.frame(OD=c(0.3,0.41))


my.dat<-read.csv("140709_diff.csv",sep="\t",header=F,row.names=letters[1:8])

measured.data<-data.frame(OD=c(my.dat[,1],my.dat[,2],my.dat[,3],my.dat[,4],my.dat[,5],my.dat[,6],my.dat[,7],my.dat[,8],my.dat[,9],my.dat[,10],my.dat[,11],my.dat[,12]))

row.names(measured.data)<-paste(rep(letters[1:8],12),rep(1:12,each=8),sep="")

#predict(FourP,measured.data)

conc<-data.frame(concentration=predict(FourP,measured.data),OD=measured.data$OD)
row.names(conc)<-paste(rep(letters[1:8],12),rep(1:12,each=8),sep="")

conc.results<-as.data.frame(matrix(conc$concentration,nrow=8),row.names=letters[1:8])
write.csv(conc.results,file="results.csv")

plot(conc$OD,conc$concentration,log="xy")
plot(FourP,add=T,col="blue",type="none")

lin.res<-read.csv("linear_results.csv",sep="\t",header=F,row.names=letters[1:8],dec=",")