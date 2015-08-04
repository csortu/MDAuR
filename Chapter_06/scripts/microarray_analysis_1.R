datapath <- '/data/work/Courses/Microarray/GSE29797'

setwd(datapath)

library(affy)
dat<-ReadAffy()

dat #this will install the first time htmg430pmcdf package
annotation(dat)
cdfName(dat)

exp.des <- pData(dat)
exp.des$Genotype <- factor(rep(c("WT","dTsc1"),each=6))
exp.des$Stimulation <- factor(rep(c("0h","4h"),6))
exp.des

pData(dat) <- exp.des

image(dat[,1])

# rudimentary QC

deg <- AffyRNAdeg(dat)
summaryAffyRNAdeg(deg)
plotAffyRNAdeg(deg)

library(affyQCReport)
correlationPlot(dat)
borderQC1(dat)



# normalization

hist(dat)
boxplot(log2(exprs(dat)))

datrma<-rma(dat)
datrma # this is an expression set

matexp<-exprs(datrma)

par(mfrow=c(1,2))
boxplot(log2(exprs(dat)),main="Before normalization",ylab="log2(Intensity)")
boxplot(matexp,main="RMA normalizaed data")
par(mfrow=c(1,1))

MAplot(datrma[,1:4],pairs=TRUE,plot.method="smoothScatter")



#datmas <- mas5(dat)
#datmas


# quinck pca for finding grouped samples

library(affycoretools)
plotPCA(matexp, groups = as.numeric(pData(dat)[,2]), groupnames = levels(pData(dat)[,2]))
plotPCA(matexp, groups = as.numeric(pData(dat)[,3]), groupnames = levels(pData(dat)[,3]))

# finding differentially expressed genes

library(limma)

genotype<- factor(pData(dat)[,2] , levels = levels(pData(dat)[,2]))
stimul <- factor(pData(dat)[,3] , levels = levels(pData(dat)[,3]))
design<- model.matrix(~genotype)


fit <-lmFit(matexp,design)
fit

fit<-eBayes(fit)

dg.top.50 <- topTable(fit, coef = 2, adjust = "fdr", n = 50)

#mouse.annot <- read.delim("GPL11180-26917.txt",skip=16,row.names=1)

mouse.annot <- read.table("GPL11180_selected_annotation.txt")
sel.annot <- mouse.annot[row.names(dg.top.50),]

dg.top.50.annotated <- cbind(dg.top.50,sel.annot)

write.table(dg.top.50.annotated[order(dg.top.50.annotated$logFC),c(1,5,7,17:15)],"diff_genes_genotype.txt")


selected <- p.adjust(fit$p.value[, 2],method="fdr") <0.03

matexp.gen <- matexp[selected,]
dim(matexp.gen)

heatmap(matexp.gen)


design.1 <-  model.matrix(~stimul)

fit.1 <-lmFit(matexp,design.1)
fit.1<-eBayes(fit.1)

dg.top.50.stim <- topTable(fit.1, coef = 2, adjust = "fdr", n = 50) 
sel.annot <- mouse.annot[row.names(dg.top.50.stim),]
dg.top.50.stim.annotated <- cbind(dg.top.50.stim,sel.annot)

write.table(dg.top.50.stim.annotated[order(dg.top.50.annotated$logFC),c(1,5,7,17:15)],"diff_genes_stimulation.txt")

dg.top.50.annotated$Gene.Symbol
dg.top.50.stim.annotated$Gene.Symbol

