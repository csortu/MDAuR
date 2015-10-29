library(mzR)
msfile<-'000_select.mzXML'


my.data<-openMSfile(msfile)
runInfo(my.data)
instrumentInfo(my.data)

header(my.data)

peaks(my.data)
peaks(my.data,1)

plot(peaks(my.data,1),type="h")

plist<-peaks(my.data,7)
plot(plist, type="h", xlab="m/z",ylab="Intensity",main="Scan 7")

close(my.data)

