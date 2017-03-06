library(GenomicAlignments)

projectDir <- getwd() # set the directory here to where you opened /home/ortutay/tmp/1
alnFiles <- list.files(path=projectDir,pattern="*bowtie2.bam",full.names=T)

## Note: if you are using Bioconductor version 14, paired with R 3.1, you should also load the following library. You do not need to load this library, and it will not be available to you, if you are using Bioconductor version 13, paired with R 3.0.x.
## library(GenomicAlignments)

aln <- readGAlignments(alnFiles)
aln

table(strand(aln))

table(width(aln))

table(cigar(aln))

head(sort(table(cigar(aln)), decreasing=TRUE))

table(seqnames(aln)) # guess which chromosome has BRCA gene...

xtabs(~seqnames + strand, as.data.frame(aln))


aln2 <- readGAlignmentPairs(alnFiles)
aln2

table(seqnames(aln2))

# see the distribution of original fragment lengths
isProperPair(aln2)
aln2.p<-aln2[isProperPair(aln2)]
plot(density(end(ranges(right(aln2.p)))-start(ranges(left(aln2.p)))),main="Fragment length distribution")


cov2 <- coverage(aln2)
cov2

mean(cov2$chr17)
max(cov2$chr17)

large.range <- slice(cov2$chr17, lower=10)
large.range.region <- data.frame(position=start(large.range):end(large.range),coverage=as.integer(large.range[[1]]))


plot(large.range.region, type="l",main="Genomic region with over 10x coverage",ylab="Coverage",sub="Chromosome 17")

hq.ranges <- slice(cov2$chr17, lower=25)

hq.longest.region <-  data.frame(position=start(hq.ranges)[7]:end(hq.ranges)[7],coverage=as.integer(hq.ranges[[7]]))
plot(hq.longest.region, type="l",main="Genomic region with over 10x coverage",ylab="Coverage",sub="Chromosome 17")

#hq.ranges.width <- width(IRanges(start=start(hq.ranges),end=end(hq.ranges)))

plot(as.integer(hq.ranges[[1]]), type="l",main="Coverage of peak 1")
plot(as.integer(hq.ranges[[10]]), type="l",main="Coverage of peak 10")




#peaks <- GRanges(seqnames=paste('peak',1:length(hq.ranges),sep="."),ranges=IRanges(start=start(hq.ranges),end=end(hq.ranges)),strand='-')
peaks <- GRanges(seqnames='chr17',ranges=IRanges(start=start(hq.ranges),end=end(hq.ranges)),strand='*')


library(TxDb.Hsapiens.UCSC.hg18.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene

seqlevels(txdb)
seqlevels(txdb, force=TRUE) <- c("chr17")


genelist <- genes(txdb)

genehits<-findOverlaps(peaks,genelist)

#genelist[unique(genehits@subjectHits)] #earlier versions used subjectHits slot for marking the covered genes
genelist[unique(genehits@to)]

gids <- genelist[unique(genehits@to)]$gene_id

#window(genelist, start=start(large.range),end=end(large.range))
cols <- c("GENEID","TXNAME","TXCHROM","TXSTRAND","TXSTART","TXEND")

gene.info<-select(txdb, gids, cols, keytype="GENEID")

vals <- list(gene_id="672")
gene.range <- genes(txdb,filter=vals)


library(ggbio)
p.ideo <- Ideogram(genome = "hg18",which=gene.range)
p.ideo + xlim(GRanges("chr17", IRanges(min(start(ranges(peaks))), max(end(ranges(peaks))))))
#gene.mods <- autoplot(txdb,which=gene.range)

#p.ideo + xlim(GRanges("chr17", IRanges(min(start(ranges(peaks))), max(end(ranges(peaks))))))

gene.mods <- autoplot(txdb,which=gene.range)

aln.s <- aln[as.character(seqnames(aln))=="chr17"] # keep only chr17 reads
read.starts <- start(aln.s)
read.ends <- end(aln.s)

aln.s.range<-aln.s[read.starts > min(start(ranges(gene.range))) & read.ends<max(end(ranges(gene.range)))]

reads <- autoplot(aln.s.range,which=gene.range)

trks <- tracks(Ideogram=p.ideo + xlim(GRanges("chr17", IRanges(min(start(ranges(peaks))), max(end(ranges(peaks)))))),Genes=gene.mods,Reads=reads,main="Reads in thse best covered region")

trks


library(VariantAnnotation)

vcfFile <- list.files(path=projectDir,pattern="*.vcf",full.names=T)
vcf.data <- readVcf(vcfFile, "hg18")
vcf.range <- as(vcf.data, "VRanges")

snps <- autoplot(vcf.range, which = gene.range, geom = "rect", arrow = FALSE)
snps
trks <- tracks(Ideogram=p.ideo + xlim(GRanges("chr17", IRanges(min(start(ranges(peaks))), max(end(ranges(peaks)))))),
               Genes=gene.mods,
               Reads=reads,
#               SNPs=snps, # plotting SNPs is broken in ggbio right now
               main="Reads in thse best covered region")

trks
