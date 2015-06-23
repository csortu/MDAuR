
a <- rep(NA,30)

# generate reference samples
a[1:3] <- rnorm(3,mean=1,sd=0.05)
for(i in 1:5){
  for(n in 1:3){
#    print(i*3+n)
#    print(i*3+n-3)
    a[i*3+n] <- a[i*3+n-3]/10
  }
}

# generate random dilutions

b <- runif(4,min=0.2,max=0.8)/10^runif(4,min=0,max=5)

for(i in 1:4){
#  print((18+i*3-2):(18+i*3))
  a[(18+i*3-2):(18+i*3)] <- b[i]*rnorm(3,mean=1,sd=0.05)
#  print(a)
}

# scale the measures to look as realistic measures
a <- a*4
# run the reactions

eff <- rnorm(length(a),mean=1.7,sd=0.02)
thres <- rnorm(length(a),mean=10000,sd=300)

dat <- matrix(nrow=40,ncol=30)

dat[1,] <- a

for(cyc in 2:40){
  for (sam in 1:length(a)){
    # implement the logistic curve
    dat[cyc,sam] <- dat[cyc-1,sam]*eff[sam]*(thres[sam]-dat[cyc-1,sam])/thres[sam]
    
  }
}

dat <- as.data.frame(dat)
names(dat) <- c(paste("RefDil",rep(0:5,each=3),".",rep(1:3,6),sep=""),paste("S",rep(1:4,each=3),".",rep(1:3,4),sep=""))


matplot(dat,t="l",col=rep(1:10,each=3),lty=c(rep(1,18),rep(2,12)))

cyc <- data.frame(Cycles=1:40)
dat <- cbind(cyc,dat)

write.table(dat,file="abs_quant_data.txt")
