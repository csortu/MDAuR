#data.path <- paste0(system.file(package = "EatonEtAlChIPseq"),"/extdata")
data.path <- dir(system.file(package = "EatonEtAlChIPseq"),pattern='extdata',full.names=T)
map.files <- list.files(data.path,pattern=".gz")
map.files

#data.path <- '/home/ortutay/tmp/HTDA/T6_data/chipseq'

library(ShortRead)

chip.aln <- readAligned(data.path,pattern="GSM424494_wt_G2_orc_chip_rep1_S",type="MAQMapview")

chip.aln

head(chromosome(chip.aln))
levels(chromosome(chip.aln))

new.chr <- chromosome(chip.aln)
levels(new.chr) <- "chrXIV"
chip.aln <- renew(chip.aln,chromosome=new.chr)
head(chromosome(chip.aln))


chip.ranges <- as(chip.aln,"GRanges")
ranges(chip.ranges)
mean(width(ranges(chip.ranges)))

library(chipseq)

estimate.mean.fraglen(chip.ranges,method ="coverage")
estimate.mean.fraglen(chip.ranges,method ="correlation")

chip.ext.ranges <- resize(chip.ranges,width=200)

chip.cov <- coverage(chip.ext.ranges)

## plot(as.numeric(dimnames(table(chip.cov))[[2]]),log2(as.numeric(table(chip.cov))),xlim=c(0,3000),xlab="Coverage",ylab="Number of bases (log2)")
## plot(as.numeric(dimnames(table(chip.cov))[[2]]),log2(as.numeric(table(chip.cov))),xlim=c(0,300),xlab="Coverage",ylab="Number of bases (log2)")

peakCutoff(chip.cov, fdr = 0.003)

chip.peaks <- slice(chip.cov,lower=8)

chip.peaks$chrXIV

peak.ranges <- peakSummary(chip.peaks)
peak.ranges

peak.ranges.df <- as.data.frame(peak.ranges)
max(peak.ranges.df$max)

rougepeak<-which(max(peak.ranges.df$max)==peak.ranges.df$max)
peak.ranges.df <- peak.ranges.df[-rougepeak,]

plot(peak.ranges.df$maxpos,peak.ranges.df$max,type="s",main="ChIP peaks\nChromosome XIV",xlab="Position",ylab="Maximum coverage",col="blue")

cov.pos <- coverage(chip.ext.ranges[strand(chip.ext.ranges)=='+'])
cov.neg <- coverage(chip.ext.ranges[strand(chip.ext.ranges)=='-'])

peaks.pos <- Views(cov.pos, ranges(chip.peaks))
peaks.neg <- Views(cov.neg, ranges(chip.peaks))

peaks.pos$chrXIV

coverageplot(peaks.pos$chrXIV[1],peaks.neg$chrXIV[1])
coverageplot(peaks.pos$chrXIV[5],peaks.neg$chrXIV[5])

## peak.ranges.df <- as.data.frame(peak.ranges)
## peak.ranges.df[which(peak.ranges.df$max>1000),]


chip.high.peaks <- slice(chip.cov,lower=150)

chip.sel.peaks <- peakSummary(chip.high.peaks$chrXIV[width(peakSummary(chip.high.peaks))>100])

bind.sites <- GRanges(seqnames='chrXIV',ranges=unlist(ranges(chip.sel.peaks)),strand='*')

# reference data from the article

peak.files <- list.files(data.path,pattern=".bed",full.names=T)
ref.peaks <- read.delim(peak.files[1],header=F,col.names=c("chr","start","end","name","p.value"))
ref.peak.ranges <- GRanges(seqnames='chrXIV',ranges=IRanges(start=ref.peaks$start,end=ref.peaks$end),strand='*')

as.data.frame(findOverlaps(bind.sites,ref.peak.ranges)) # great! it is almost identical

# but which genes are there?

library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene

seqlevels(txdb)
seqlevels(txdb, force=TRUE) <- c("chrXIV")
genelist <- genes(txdb)

covers <- findOverlaps(bind.sites,genelist,ignore.strand=T)

affected.genes <- as.data.frame(genelist[subjectHits(covers)])

library(org.Sc.sgd.db)

select(org.Sc.sgd.db, as.data.frame(affected.genes)$gene_id, keytype="ORF",columns=c("GENENAME","ENTREZID"))

