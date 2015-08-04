library(GEOquery)

datadir <- '/home/data/work/Courses/Microarray/GSE43805/GEO'

gse <- getGEO("GSE43805",destdir=datadir)
gse
gse[[1]]
illumina.dat.geo <- gse[[1]]


pData(phenoData(illumina.dat.geo))[c(1,2,8)] #34 columns

treat <- pData(phenoData(illumina.dat.geo))[,8]
design <- model.matrix(~treat)

library(limma)

fit.geo <-lmFit(illumina.dat.geo,design)
fit.geo<-eBayes(fit.geo)

dg.top.50.geo <- topTable(fit.geo, coef = 1, adjust = "fdr", n = 50)
write.table(dg.top.50.geo[,c(8,10:13,30:36)],file="illumina_annotated_dg.txt")

selected <- p.adjust(fit.geo$p.value[, 2],method="fdr") <0.03
sign.genes <- exprs(illumina.dat.geo)[selected,]

heatmap(sign.genes,labCol=c("Control","Control","Virus691","Virus94"),margins=c(10,5))

library(gplots) 
my.palette <-rev(redgreen(75))

heatmap.2(sign.genes,labCol=c("Control","Control","Virus691","Virus94"),margins=c(10,5),col=my.palette)

