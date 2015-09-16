library(RNAseqData.HNRNPC.bam.chr14)
alnFiles <- RNAseqData.HNRNPC.bam.chr14_BAMFILES

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

seqlevels(txdb, force=TRUE) <- c("chr14")
genelist <- genes(txdb)

library(GenomicAlignments)
counts <- data.frame(row.names=genelist$gene_id)

for (i in 1:length(RNAseqData.HNRNPC.bam.chr14_BAMFILES)){
  print(paste("Reading:",alnFiles[i]))
  aln <- readGAlignments(alnFiles[i])
  reads <- GRanges(seqnames='chr14',ranges=ranges(aln),strand='*')
  genecounts<-countOverlaps(genelist, reads)
  a <- data.frame(genecounts)
  names(a) <- RNAseqData.HNRNPC.bam.chr14_RUNNAMES[i]
  counts <- cbind(counts,a)
}

counts <- counts[!rowSums(counts)==0,]

library(org.Hs.eg.db)
sym <- select(org.Hs.eg.db, row.names(counts), keytype="ENTREZID","SYMBOL")
row.names(sym) <- sym$ENTREZID
sym$SYMBOL[is.na(sym$SYMBOL)]<-sym$ENTREZID[is.na(sym$SYMBOL)] # If there are no matching gene symbol available, gene ID is used.
counts <- cbind(sym$SYMBOL,counts)
row.names(counts) <- as.character(counts$"sym$SYMBOL")
counts$"sym$SYMBOL" <- NULL

## table 7.1

write.csv(counts[c("HNRNPC","HSPA2"),],file="table7_1.csv")

## Figure 7.1

counts <- read.table("HNRNPC_chr14_gene_read_counts.txt")
exp.des <- factor(rep(c('Control','HNRNPC knockdown'),each=4))
names(exp.des) <- RNAseqData.HNRNPC.bam.chr14_RUNNAMES

library(edgeR)
genexp <- DGEList(counts=counts,group=exp.des)
genexp <- calcNormFactors(genexp)
genexp <- estimateCommonDisp(genexp)
genexp <- estimateTagwiseDisp(genexp)

genediff <- exactTest(genexp)
genediff$table <- cbind(genediff$table,FDR=p.adjust(genediff$table$PValue,method="fdr"))
gsign <- genediff$table[genediff$table$FDR<0.03,]
gsign <- gsign[order(gsign$FDR),]

write.csv(format(head(gsign),digits=2),file="table7_2.csv")


library(DESeq)
count.set <- newCountDataSet(counts,exp.des)

count.set <- estimateSizeFactors(count.set)

count.set <- estimateDispersions(count.set)

tiff("figure7_1.tif",width=5,height=5,units="in",res=1200)
plot(cbind(genexp$tagwise.dispersion,fData(count.set)),
     ylim=c(0,2),
     main="Gene-wide dispersion",xlab="edgeR",ylab="DESeq",
     pch=20)
abline(0,1,lty=2)
dev.off()

## figure 7.2

genediff.ds <- nbinomTest(count.set,"Control","HNRNPC knockdown")

#plotMA(genediff.ds)
tiff("figure7_2.tif",width=5,height=5,units="in",res=1200)
plot(log2FoldChange ~ baseMean,data=genediff.ds[genediff.ds$padj>0.1,],log="x",
     xlab="Count number",
     ylab="log fold change",
     main="MA plot",
     pch=20)
abline(h=0,lwd=2)
points(genediff.ds[genediff.ds$padj<0.1,c("baseMean","log2FoldChange")],
       pch=8)
legend(1800,3,pch=8,legend = "significant\nDE",bty="n")
dev.off()

## figure 7.3
tiff("figure7_3.tif",width=5,height=5,units="in",res=1200)
par(mfrow=c(1,3),mai=c(0.52,0.67,0.42,0.12))
plot(gsign[intersect(row.names(gsign),gsign.ds$id),"logFC"],gsign.ds[intersect(row.names(gsign),gsign.ds$id),"log2FoldChange"],main="Log Fold Changes",xlab="edgeR",ylab="DESeq")
plot(gsign[intersect(row.names(gsign),gsign.ds$id),"FDR"],gsign.ds[intersect(row.names(gsign),gsign.ds$id),"padj"],main="Adj. p-values",xlab="edgeR",ylab="DESeq")
plot(gsign[intersect(row.names(gsign),gsign.ds$id),"FDR"],gsign.ds[intersect(row.names(gsign),gsign.ds$id),"padj"],main="Adj. p-values",xlab="edgeR",ylab="DESeq",log="xy",xlim=c(1e-20,1e-1),ylim=c(1e-20,1e-1))
dev.off()

## figure 7.4

#par(mfrow=c(1,1))

library(gplots) 
my.palette <-rev(redgreen(75))

tiff("figure7_4_color.tif",width=5,height=5,units="in",res=1200)
heatmap.2(log2(as.matrix(counts[rownames(gsign),]+1)),
          col=my.palette,
          margins=c(8,1),
          key = F)
dev.off()
