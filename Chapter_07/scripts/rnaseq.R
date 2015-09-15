# Get the data

library(RNAseqData.HNRNPC.bam.chr14)

exp.des <- factor(rep(c('Control','HNRNPC knockdown'),each=4))
names(exp.des) <- RNAseqData.HNRNPC.bam.chr14_RUNNAMES



alnFiles <- RNAseqData.HNRNPC.bam.chr14_BAMFILES


# Annotation

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb)
seqlevels(txdb, force=TRUE) <- c("chr14")
genelist <- genes(txdb)
txlist <- transcripts(txdb)
exonlist <- exons(txdb)

genelist
txlist
exonlist

# small game of alternative transcript statistics

alt.trans<-findOverlaps(txlist,genelist)
plot(table(table(alt.trans@subjectHits)))
table(table(alt.trans@subjectHits))

# Calculate counts per gene

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
head(counts)
dim(counts)

# eliminate genes without any reads

sum(rowSums(counts)==0)
sum(rowSums(counts)<10)

counts <- counts[!rowSums(counts)==0,]
#counts1 <- counts
# Add Symbols as row names

library(org.Hs.eg.db)

sym <- select(org.Hs.eg.db, row.names(counts), keytype="ENTREZID","SYMBOL")
row.names(sym) <- sym$ENTREZID
sym$SYMBOL[is.na(sym$SYMBOL)]<-sym$ENTREZID[is.na(sym$SYMBOL)]
counts <- cbind(sym$SYMBOL,counts)
row.names(counts) <- as.character(counts$"sym$SYMBOL")
counts$"sym$SYMBOL" <- NULL

# it is wise to save at this point the table as a text file
counts["HNRNPC",]
counts["HSPA2",]

#t.test(counts["SPTLC2",1:4],counts["SPTLC2",5:8],alternative="less")

write.table(counts,file="HNRNPC_chr14_gene_read_counts.txt")

stop("data is ready")
# Differential expression calculation with edgeR

counts <- read.table("HNRNPC_chr14_gene_read_counts.txt")

# naive t test approach

counts["SPTLC2",]

t.test(counts["SPTLC2",1:4],counts["SPTLC2",5:8])

tt.p.vals <- apply(counts,1,function(dat){t.test(x = dat[1:4], y = dat[5:8])$p.value})
tt.adj.p.vals <- p.adjust(tt.p.vals,method="fdr")
counts <- cbind(counts,tt.p.vals,tt.adj.p.vals)


head(counts[order(tt.adj.p.vals),])

counts <- counts[order(tt.adj.p.vals),]
gsign.tt <- counts[counts$tt.adj.p.vals<0.03,]

dim(gsign.tt)
gsign.tt["BDKRB1",]


# exact test with edgeR

counts <- read.table("HNRNPC_chr14_gene_read_counts.txt")


library(edgeR)

# we should filter genes with lower than 5 

# classic

genexp <- DGEList(counts=counts,group=exp.des)
genexp <- calcNormFactors(genexp)
genexp <- estimateCommonDisp(genexp)
genexp <- estimateTagwiseDisp(genexp)

#head(genexp$counts[order(genexp$tagwise.dispersion),])
#tail(genexp$counts[order(genexp$tagwise.dispersion),])

genediff <- exactTest(genexp)
topTags(genediff)

genediff$table <- cbind(genediff$table,FDR=p.adjust(genediff$table$PValue,method="fdr"))

gsign <- genediff$table[genediff$table$FDR<0.03,]
gsign <- gsign[order(gsign$FDR),]
dim(gsign)
head(gsign)

# general linear models

repls <- rep(rep(1:2,each=2),2)
design <- model.matrix(~exp.des+repls)

genexp.glm <- DGEList(counts=counts)
genexp.glm <- estimateGLMCommonDisp(genexp.glm,design)
#genexp.glm <- estimateGLMTrendedDisp(genexp.glm,design)
genexp.glm <- estimateGLMTagwiseDisp(genexp.glm,design)
gene.exp.fit <- glmFit(genexp.glm,design)
genediff.lrt <- glmLRT(gene.exp.fit,coef=2)
topTags(genediff.lrt)

genediff.lrt.2 <- glmLRT(gene.exp.fit,coef=3)
topTags(genediff.lrt.2)

genediff.lrt$table <- cbind(genediff.lrt$table,FDR=p.adjust(genediff.lrt$table$PValue,method="fdr"))

gsign.lrt <- genediff.lrt$table[genediff.lrt$table$FDR<0.03,]
gsign.lrt <- gsign.lrt[order(gsign.lrt$FDR),]


dim(gsign.lrt)
head(gsign.lrt)
head(gsign)

# differential gene expression calculation with DESeq

library(DESeq)

count.set <- newCountDataSet(counts,exp.des)


count.set <- estimateSizeFactors(count.set)
sizeFactors(count.set)

plot(genexp$samples$norm.factors,sizeFactors(count.set),xlab="edgeR",ylab="DESeq",main="Normalization factors",type="n")
abline(0,1,lty=2)
text(genexp$samples$norm.factors,sizeFactors(count.set),labels=names(sizeFactors(count.set)),cex=0.8)

# dispersion estimation

count.set <- estimateDispersions(count.set)
par(mfrow=c(1,3))
plotDispEsts(count.set,main="DESeq: genewise dispersion",cex=1)
plotBCV(genexp,main="edgeR: genewise dispersion normalized counts",cex=1)
plotBCV(genexp.glm,main="edgeR: genewise dispersion trend for GLM",cex=1)
par(mfrow=c(1,1))

# The dispersion values for individual genes
head(cbind(genexp$tagwise.dispersion,fData(count.set)))
plot(cbind(genexp$tagwise.dispersion,fData(count.set)),ylim=c(0,2),main="Gene-wise dispersion",xlab="edgeR",ylab="DESeq")
abline(0,1,lty=2)


# simple differential expression

genediff.ds <- nbinomTest(count.set,"Control","HNRNPC knockdown")
plotMA(genediff.ds)

gsign.ds <- genediff.ds[genediff.ds$padj<0.03,]
gsign.ds <- gsign.ds[order(gsign.ds$padj),]
row.names(gsign.ds)<-gsign.ds$id

head(gsign.ds,10)
head(gsign,10)

par(mfrow=c(1,3))
plot(gsign[intersect(row.names(gsign),gsign.ds$id),"logFC"],gsign.ds[intersect(row.names(gsign),gsign.ds$id),"log2FoldChange"],main="Log Fold Changes",xlab="edgeR",ylab="DESeq")
plot(gsign[intersect(row.names(gsign),gsign.ds$id),"FDR"],gsign.ds[intersect(row.names(gsign),gsign.ds$id),"padj"],main="Adjusted p-values",xlab="edgeR",ylab="DESeq")
plot(gsign[intersect(row.names(gsign),gsign.ds$id),"FDR"],gsign.ds[intersect(row.names(gsign),gsign.ds$id),"padj"],main="Adjusted p-values",xlab="edgeR",ylab="DESeq",log="xy",xlim=c(1e-20,1e-1),ylim=c(1e-20,1e-1))
par(mfrow=c(1,1))

# heatmap visualization

diff.exp <- data.frame(edgeR=gsign[intersect(row.names(gsign),gsign.ds$id),"logFC"],DESeq=gsign.ds[intersect(row.names(gsign),gsign.ds$id),"log2FoldChange"])
row.names(diff.exp) <- intersect(row.names(gsign),gsign.ds$id)

diff.exp[is.infinite(diff.exp$DESeq),]<--2

library(gplots) 
my.palette <-rev(redgreen(75))

# full dataset
heatmap.2(log2(as.matrix(counts+1)),col=my.palette,margins=c(10,12))

# Behaviour of differentially expressed genes
heatmap.2(as.matrix(diff.exp),Colv=F,col=my.palette,margins = c(12, 10))


# multifactor fifferential expression with glm

exp.des.cmplx<-data.frame(Type=as.character(exp.des),Repl=repls)
row.names(exp.des.cmplx)<-names(exp.des)

count.set.glm <- newCountDataSet(counts,exp.des.cmplx)

count.set.glm <- estimateSizeFactors(count.set.glm)
count.set.glm <- estimateDispersions(count.set.glm)

# fit models to count data

genediff.ds.mod1 <- fitNbinomGLMs(count.set.glm, count ~ Type + Repl)
genediff.ds.mod2 <- fitNbinomGLMs(count.set.glm, count ~ Type)

genediff.ds.glm.pvals<-nbinomGLMTest(genediff.ds.mod1,genediff.ds.mod2)
genediff.ds.glm.fdr<-p.adjust(genediff.ds.glm.pvals,method="fdr")

genediff.ds.glm <- cbind(genediff.ds.mod1[,2:3],PVal=genediff.ds.glm.pvals,FDR=genediff.ds.glm.fdr)

gsign.ds.glm <- genediff.ds.glm[genediff.ds.glm$FDR<0.03,]
gsign.ds.glm <- gsign.ds.glm[order(gsign.ds.glm$FDR),]

# compare to edgeR

common.genes <- intersect(row.names(gsign.lrt),row.names(gsign.ds.glm))
dim(gsign.lrt)
dim(gsign.ds.glm)

plot(gsign.lrt[common.genes,"logFC"],gsign.ds.glm[common.genes,"TypeHNRNPC knockdown"],main="Fold cahnge ith GLMs",xlab="edgeR",ylab="DESeq")
abline(0,1,lty=2)



