library(Cairo)

my.phred <- function(x){-10*log10(x)}

CairoTIFF("figure4_2.tif",width=5,height=5,units="in",res=1200)

curve(my.phred,0.001,1,main="Sanger quality score",xlab="Probability of incorrect base",ylab="Quality score")
text(0.22,15,labels="Q=-10 x log(p)")
abline(v=0.05,lwd=0.5)
text(0.11,0,labels="0.05")

dev.off()



