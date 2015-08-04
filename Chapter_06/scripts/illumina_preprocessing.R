datpath <- '/data/work/Courses/Microarray/GSE43805/'

setwd(datpath)

dat <- read.delim("GSE43805_non_normalized.txt",row.names=1)

exprs <- as.matrix(dat[,c(1,3,5,7)])

library(limma)

exprs.norm <- normalizeBetweenArrays(log2(exprs))
exprs.norm.loess <- normalizeBetweenArrays(log2(exprs),method="cyclicloess")

library(Biobase)
minimalSet <- ExpressionSet(assayData=exprs.norm)

pheno <- data.frame(treatment=c("control","control","virus691","virus94"),row.names=colnames(exprs))

metadata <- data.frame(labelDescription="Lentivirus treatment of samples",row.names="treatment")

phenoData <- new("AnnotatedDataFrame",data=pheno, varMetadata=metadata)

illumina.dat <- ExpressionSet(assayData=exprs.norm,phenoData=phenoData,annotation="Illumina HumanHT-12 V4.0 expression beadchip")

treat <- pData(illumina.dat)[,1]
design <- model.matrix(~treat)


fit <-lmFit(illumina.dat,design)
fit<-eBayes(fit)

dg.top.50.v691 <- topTable(fit, coef = 1, adjust = "fdr", n = 50) 
dg.top.50.v94  <- topTable(fit, coef = 2, adjust = "fdr", n = 50) 

