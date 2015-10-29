library(mzR)
msfile<-'000_select.mzXML'

my.data<-openMSfile(msfile)

pl<-peaks(my.data,7)

#pl.mz<-pl[,1]
#pl15<-pl[pl[,2]%in% head(sort(pl[,2],decreasing=T),15),1]

topnum<-150
prec<-0.04
pl.top<-pl[pl[,2]%in% head(sort(pl[,2],decreasing=T),topnum),1]

peakdiff<-outer(pl.top,pl.top,'-')
#plot(density(abs(peakdiff)))

#plot(pl,type="h",xlim=c(300,600),ylim=c(0,500000))

peakdiff.aa<-matrix(rep('-',topnum^2),ncol=topnum,nrow=topnum)

aa.mass<-read.csv("aa_mass.csv",row.names=1)

ystep<-max(pl[,2])/25
ypos<-ystep*4
plot(pl,type="h")
for (aa in row.names(aa.mass)){

  peakdiff.aa[abs(peakdiff-aa.mass[aa,])<prec]<-aa              
 
  aa.match<-which(abs(peakdiff-aa.mass[aa,])<prec,arr.ind=T)
  if(length(aa.match)){
    for (i in 1:dim(aa.match)[1]){
      segments(pl.top[aa.match[i,2]],ypos,pl.top[aa.match[i,1]],ypos)
#      print(pl.top[aa.match[i,]])
       text(pl.top[aa.match[i,2]],ypos,aa)

    }
  }

  ypos<-ypos+ystep

}
