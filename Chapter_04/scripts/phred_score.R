my.phred <- function(x){-10*log10(x)}
curve(my.phred,0.001,1,main="Sanger quality score",xlab="Probability of incorrect base",ylab="Quality score")
text(0.2,15,labels="Q=-10 x log(p)")
abline(v=0.05,lty=2)
text(0.1,0,labels="0.05")

my.odd <- function(x){x/(1-x)}
my.solexa <- function(x){my.phred(my.odd(x))}
curve(my.solexa,0,1,add=T,col="red")