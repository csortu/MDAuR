library(TxDb.Hsapiens.UCSC.hg18.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene
seqlevels(txdb, force=TRUE) <- c("chr17")


library(GenomicAlignments)

projectDir <- '/home/ortutay/tmp/HTDA/Topic2'
alnFiles <- list.files(path=projectDir, pattern="*.bam", full.names=T)
aln <- readGAlignments(alnFiles)

aln2 <- readGAlignmentPairs(alnFiles)
cov2 <- coverage(aln2)
hq.ranges <- slice(cov2$chr17, lower=25)
peaks <- GRanges(seqnames='chr17', ranges=IRanges(start=start(hq.ranges), end=end(hq.ranges)), strand='*')




library(ggbio)
vals <- list(gene_id="672")
gene.range <- genes(txdb,vals)
p.ideo <- Ideogram(genome = "hg18",which=gene.range)

gene.mods <- autoplot(txdb,which=gene.range,names.expr="",gap.geom = "chevron")

aln.s <- aln[as.character(seqnames(aln))=="chr17"]
read.starts <- start(aln.s)
read.ends <- end(aln.s)
aln.s.range<-aln.s[read.starts > min(start(ranges(gene.range))) & read.ends<max(end(ranges(gene.range)))]
reads <- autoplot(aln.s.range,which=gene.range)



trks <- tracks(Ideogram=p.ideo + xlim(GRanges("chr17", IRanges(min(start(ranges(peaks))), max(end(ranges(peaks)))))),Transcripts=gene.mods,Reads=reads,main="Reads in the best covered region")

library(Cairo)

CairoTIFF("figure4_4col.tif",width=5,height=5,units="in",res=1200)

trks

dev.off()
