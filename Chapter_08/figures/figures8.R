data.dir <- '/home/data/work/Courses/Microarray/ChIP-chip/E-GEOD-31301'
setwd(data.dir)
my.dat<-read.delim("GSM775491_sample_table.txt",row.names=1)
adf <- read.delim('A-GEOD-13972.adf.txt',skip=15,row.names=1)
names(adf) <- c("Coordinate","Sequence")
adf$chr <- apply(matrix(adf$Coordinate,ncol=1),1,function(coord){strsplit(coord,':')[[1]][1]})
adf$pos <- apply(matrix(adf$Coordinate,ncol=1),1,function(coord){strsplit(coord,':')[[1]][2]})
adf$start <- apply(matrix(adf$pos,ncol=1),1,function(coord){strsplit(coord,'-')[[1]][1]})
adf$end <- apply(matrix(adf$pos,ncol=1),1,function(coord){strsplit(coord,'-')[[1]][2]})
adf$start <- as.numeric(adf$start)
adf$end <- as.numeric(adf$end)
adf$chr <- as.factor(adf$chr)
adf <- adf[!is.na(adf$chr),]
my.dat <- cbind(my.dat,adf[,c(3,5,6)])

## figure 8.1

#plot(density(my.dat$VALUE,na.rm=T))
tiff("figure8_1.tif",width=5,height=5,units="in",res=1200)
plot(density(my.dat$VALUE,na.rm=T),
     xlim=c(-2,2),xlab="Ratio",
     main="Distribution of immunoprecipitation ratios")
dev.off()

## figure 8.2
my.dat.c4 <- my.dat[my.dat$chr=='chr4',]
my.dat.c4 <- my.dat.c4[order(my.dat.c4$start),]
w.start <- 500000
w.end <- 512000

tiff("figure8_2.tif",width=5,height=5,units="in",res=1200)
plot(my.dat.c4$VALUE ~ my.dat.c4$start,
     t='l',xlim=c(w.start,w.end),
     main="ChIP-Chip peaks\non a segment of Chromosome 4",
     xlab="Position",ylab="log2 IP/Input")
abline(h=0,lty=2)
segments(501000,-2.5,501300,lwd=2)
segments(507800,1.7,509800,lwd=2)
dev.off()

## figure 8.3

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

tiff("figure8_3.tif",width=5,height=5,units="in",res=1200)
plot(my.dat.c4$VALUE ~ my.dat.c4$start,
     t='l',xlim=c(w.start,w.end),
     main="ChIP-Chip peaks\non a segment of Chromosome 4"
     ,xlab="Position",ylab="log2 IP/Input")
abline(h=0,lty=2)
lines(w.pos,w.smooth,lty=3, lwd=3)
legend("topleft",bty="n",
       legend=c("raw values","smoothed values","binding region"),
       lty=c(1,3,5),lwd=c(1,3,2))
text(505507,-2.25,labels="eliminated peaks")
arrows(c(503000,503950,506926,507315),
       c(-2.247,-1.76,-1.96,-2.01),
       c(501446,503560,507426,510653),
       c(-2.23,-1.07,-1.42,-0.58),
       length=0.05, code=2)
segments(507800,1.7,509800,lty=5,lwd=2)
dev.off()

## figure 8.4

reporters <- row.names(my.dat.c4)
my.dat.2<-read.delim("GSM775492_sample_table.txt",row.names=1)
my.dat.3<-read.delim("GSM775493_sample_table.txt",row.names=1)
my.dat.c4$rep2 <- my.dat.2[reporters,]
my.dat.c4$rep3 <- my.dat.3[reporters,]

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y,use="pairwise.complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

tiff("figure8_4.tif",width=5,height=5,units="in",res=1200)
pairs(my.dat.c4[,c(1,5,6)],
      upper.panel=panel.smooth,
      lower.panel=panel.cor,
      main="Correlation of biological replicates",
      col.smooth="gray")
dev.off()

## figure 8.5

tiff("figure8_5.tif",width=5,height=5,units="in",res=1200)
plot(my.dat.c4$VALUE ~ my.dat.c4$start,t='l',xlim=c(w.start,w.end),
     main="ChIP-Chip peaks\non a segment of Chromosome 4",
     xlab="Position",ylab="log2 IP/Input")
abline(h=0,lty=2)
lines(my.dat.c4$start,my.dat.c4$rep2,lty=3,lwd=2)
lines(my.dat.c4$start,my.dat.c4$rep3,lty=5)
legend("topleft",bty="n",
       legend=c("original measurement","replicate 2","replicate 3"),
       lty=c(1,3,5),lwd=c(1,2,1))
dev.off()

## table 8.1

data.path <- dir(system.file(package = "EatonEtAlChIPseq"),pattern='extdata',full.names=T)
map.files <- list.files(data.path,pattern=".gz")
library(ShortRead)
chip.aln <- readAligned(data.path, pattern="GSM424494_wt_G2_orc_chip_rep1_S", type="MAQMapview")
new.chr <- chromosome(chip.aln)
levels(new.chr) <- "chrXIV"
chip.aln <- renew(chip.aln,chromosome=new.chr)
chip.ranges <- as(chip.aln,"GRanges")
library(chipseq)
chip.ext.ranges <- resize(chip.ranges,width=200)
chip.cov <- coverage(chip.ext.ranges)
chip.peaks <- slice(chip.cov,lower=8)
peak.ranges <- peakSummary(chip.peaks)
peak.ranges.df <- as.data.frame(peak.ranges)
write.csv(head(peak.ranges.df[order(peak.ranges.df$max,decreasing=T),]),
          file = "table8_1.csv")

## figure 8.6

cov.pos <- coverage(chip.ext.ranges[strand(chip.ext.ranges)=='+'])
cov.neg <- coverage(chip.ext.ranges[strand(chip.ext.ranges)=='-'])
peaks.pos <- Views(cov.pos, ranges(chip.peaks))
peaks.neg <- Views(cov.neg, ranges(chip.peaks))

tiff("figure8_6.tif",width=5,height=5,units="in",res=1200)
coverageplot(peaks.pos$chrXIV[47],peaks.neg$chrXIV[47],
             main="Coverage plot of a peak\non Chromosome XIV")
dev.off()

## figure 8.7

bed.file <- 'GSM1504361_PU.1_ChIPseq_PU.1_B_cells.bed.gz'
peak.tab <- read.delim(bed.file,header=F)
names(peak.tab) <- c('chr','start','end','peakName','score')

peak.tab.sel <- peak.tab[peak.tab$score %in% tail(sort(peak.tab$score),50),]


library(rGADEM)
peak.ranges.sel <- RangedData(IRanges(start=peak.tab.sel$start,
                                      end=peak.tab.sel$end),
                              space=paste0('chr',peak.tab.sel$chr))
library(BSgenome.Mmusculus.UCSC.mm9)
gadem.sel <- GADEM(peak.ranges.sel,verbose=1,genome=Mmusculus)
bestmotif <- which(nOccurrences(gadem.sel)==max(nOccurrences(gadem.sel)))
motifs.sel <- getPWM(gadem.sel)

tiff("figure8_7_color.tif",width=5,height=5,units="in",res=1200)
seqLogo(motifs.sel[[bestmotif]])
dev.off()
