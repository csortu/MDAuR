data.dir <- '/home/data/work/Courses/Microarray/ChIP-chip/E-GEOD-31301'

setwd(data.dir)

my.dat<-read.delim("GSM775491_sample_table.txt",row.names=1)

adf <- read.delim('A-GEOD-13972.adf.txt',skip=15,row.names=1)
names(adf) <- c("Coordinate","Sequence")
head(my.dat)
head(adf)

adf$chr <- apply(matrix(adf$Coordinate,ncol=1),1,function(coord){strsplit(coord,':')[[1]][1]})
adf$pos <- apply(matrix(adf$Coordinate,ncol=1),1,function(coord){strsplit(coord,':')[[1]][2]})
adf$start <- apply(matrix(adf$pos,ncol=1),1,function(coord){strsplit(coord,'-')[[1]][1]})
adf$end <- apply(matrix(adf$pos,ncol=1),1,function(coord){strsplit(coord,'-')[[1]][2]})
adf <- adf[!is.na(adf$chr),]

adf$start <- as.numeric(adf$start)
adf$end <- as.numeric(adf$end)
adf$chr <- as.factor(adf$chr)

my.dat <- cbind(my.dat,adf[,c(3,5,6)])
head(my.dat)

table(my.dat$end-my.dat$start+1)

table(my.dat$chr)

plot(density(my.dat$VALUE,na.rm=T))

my.dat.c4 <- my.dat[my.dat$chr=='chr4',]
my.dat.c4 <- my.dat.c4[order(my.dat.c4$start),]

w.start <- min(my.dat.c4$start)
w.end <- max(my.dat.c4$end)

plot(my.dat.c4$VALUE ~ my.dat.c4$start,t='l',xlim=c(w.start,w.end),main="ChIP-Chip peaks\non Chromosome 4",xlab="Position",ylab="log2 IP/Input")
abline(h=0,col="blue",lty=3)

w.start <- 500000
w.end <- 512000

plot(my.dat.c4$VALUE ~ my.dat.c4$start,t='l',xlim=c(w.start,w.end),main="ChIP-Chip peaks\non Chromosome 4",xlab="Position",ylab="log2 IP/Input")
abline(h=0,col="blue",lty=3)

w.size <- 1000

w.pos <- rep(NA,(w.end-w.start)/w.size)
w.smooth <- w.pos

i <- 1
  
for(pos in seq(w.start,w.end+w.size,w.size)){
  vals <- my.dat.c4$VALUE[my.dat.c4$start>pos &  my.dat.c4$start<(pos+w.size)]
  w.pos[i] <- pos
  w.smooth[i] <- mean(vals)
  i <- i+1
}

plot(my.dat.c4$VALUE ~ my.dat.c4$start,t='l',xlim=c(w.start,w.end),main="ChIP-Chip peaks\non Chromosome 4",xlab="Position",ylab="log2 IP/Input")
abline(h=0,col="blue",lty=3)
lines(w.pos,w.smooth,col="red",lty=2)

## replicates

reporters <- row.names(my.dat.c4)

my.dat.2<-read.delim("GSM775492_sample_table.txt",row.names=1)
my.dat.3<-read.delim("GSM775493_sample_table.txt",row.names=1)

my.dat.c4$rep2 <- my.dat.2[reporters,]
my.dat.c4$rep3 <- my.dat.3[reporters,]

plot(my.dat.c4$VALUE,my.dat.c4$rep2)
#plot(my.dat.c4$VALUE,my.dat.c4$rep2,my.dat.c4$rep3)

pairs(my.dat.c4[,c(1,5,6)])

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y,use="pairwise.complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
pairs(my.dat.c4[,c(1,5,6)],upper.panel=panel.smooth,lower.panel=panel.cor)

plot(my.dat.c4$VALUE ~ my.dat.c4$start,t='l',xlim=c(w.start,w.end),main="ChIP-Chip peaks\non Chromosome 4",xlab="Position",ylab="log2 IP/Input")
abline(h=0,col="blue",lty=3)
lines(my.dat.c4$start,my.dat.c4$rep2,col="red")
lines(my.dat.c4$start,my.dat.c4$rep3,col="green")

