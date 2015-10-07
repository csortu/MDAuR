###################################################################
# Co-expression networks

## figure 9.1

library(qpgraph)
furin.data<-read.csv("furin_significant_genes.csv",row.names=1)
pcc <- qpPCC(furin.data)
pcc$Rsign<-pcc$R
pcc$Rsign[pcc$P>0.05]<-NA


nrr.q1 <- qpNrr(as.matrix(furin.data),q=1)
nrr.q3 <- qpNrr(as.matrix(furin.data),q=3)
nrr.q5 <- qpNrr(as.matrix(furin.data),q=5)
nrr.q7 <- qpNrr(as.matrix(furin.data),q=7)

tiff("figure9_1.tif",width=5,height=5,units="in",res=1200)
par(mfrow=c(3,2),mai=c(0.42,0.42,0.42,0.12)) 
plot(density(as.numeric(pcc$R),na.rm=T),
     main="Pearson correlation\nAll",xlab="",ylab="")
mtext("A", side = 3, line = 0, adj = 0, cex = 1)
plot(density(as.numeric(pcc$Rsign),na.rm=T,xlab=""),
     main="Pearson correlation\nSignificants",xlab="",ylab="")
mtext("B", side = 3, line = 0, adj = 0, cex = 1)
plot(density(as.numeric(nrr.q1),na.rm=T),
     main="Non-rejection rate\nq=1",xlab="",ylab="")
mtext("C", side = 3, line = 0, adj = 0, cex = 1)
plot(density(as.numeric(nrr.q3),na.rm=T),
     main="Non-rejection rate\nq=3",xlab="",ylab="")
mtext("D", side = 3, line = 0, adj = 0, cex = 1)
plot(density(as.numeric(nrr.q5),na.rm=T),
     main="Non-rejection rate\nq=5",xlab="",ylab="")
mtext("E", side = 3, line = 0, adj = 0, cex = 1)
plot(density(as.numeric(nrr.q7),na.rm=T),
     main="Non-rejection rate\nq=7",xlab="",ylab="")
mtext("F", side = 3, line = 0, adj = 0, cex = 1)
dev.off()

## figure 9.2

tiff("figure9_2.tif",width=5,height=5,units="in",res=1200)
par(mfrow=c(2,2),mai=c(0.52,0.42,0.52,0.12))
plot(density(as.numeric(abs(pcc$Rsign)),na.rm=T),
     main="Pearson correlation\nSignificants",ylim=c(0,6))
plot(density(as.numeric(nrr.q1),na.rm=T),
     main="Non-rejection rate\nq=1",xlim=c(1,0),ylim=c(0,6))
plot(density(as.numeric(nrr.q3),na.rm=T),
     main="Non-rejection rate\nq=3",xlim=c(1,0),ylim=c(0,6))
plot(density(as.numeric(nrr.q5),na.rm=T),
     main="Non-rejection rate\nq=5",xlim=c(1,0),ylim=c(0,6))
dev.off()

## figure 9.3

tiff("figure9_3.tif",width=5,height=5,units="in",res=1200)
par(mfrow=c(2,2),mai=c(0.52,0.42,0.52,0.12))
qpGraphDensity(nrr.q1, title="q=1", breaks=10)
qpGraphDensity(nrr.q3, title="q=3", breaks=10)
qpGraphDensity(nrr.q5, title="q=5", breaks=10)
qpGraphDensity(nrr.q7, title="q=7", breaks=10)
dev.off()

## figure 9.4
gm <- as.matrix(nrr.q3)
thres <- 0.1
my.nodes <- row.names(gm)
edL <- vector("list", length=length(my.nodes))
names(edL) <- my.nodes
for(i in 1:length(my.nodes)){
  edL[[i]] <- list(edges=names(which(gm[i,]<thres)), weights=gm[i,which(gm[i,]<thres)])
}
library(graph)
g <- graphNEL(nodes=my.nodes, edgeL=edL)




# deprecated code
# pcc100 <- qpGraph(as(1-abs(pcc$Rsign),"dspMatrix"),
#                   topPairs=100,
#                   return.type="graphNEL")


gmpc100 <- as.matrix(abs(pcc$Rsign))
gmpc100[gmpc100==1] <- NA
topnum <- 200

pc100vals <- as.vector(gmpc100)
pc100vals <- sort(pc100vals[!is.na(pc100vals)],decreasing = T)

thres <- pc100vals[topnum]

my.nodes.pc100 <- row.names(gmpc100)
edLpc100 <- vector("list", length=length(my.nodes.pc100))
names(edLpc100) <- my.nodes.pc100
for(i in 1:length(my.nodes.pc100)){
  edLpc100[[i]] <- list(edges=names(which(gmpc100[i,]>=thres)), 
                        weights=gmpc100[i,which(gmpc100[i,]>=thres)])
}

pcc100 <- graphNEL(nodes=my.nodes.pc100, edgeL=edLpc100)

# deprecated code:
# qpg100 <- qpGraph(nrr.q3, 
#                   topPairs=100, 
#                   return.type="graphNEL")

gmpg100 <- as.matrix(nrr.q3)
gmpg100[gmpg100==1] <- NA
topnum <- 200

pg100vals <- as.vector(gmpg100)
pg100vals <- sort(pg100vals[!is.na(pg100vals)])

thres <- pg100vals[topnum]

my.nodes.pg100 <- row.names(gmpg100)
edLpg100 <- vector("list", length=length(my.nodes.pg100))
names(edLpg100) <- my.nodes.pg100
for(i in 1:length(my.nodes.pg100)){
  edLpg100[[i]] <- list(edges=names(which(gmpg100[i,]<=thres)), 
                        weights=gmpg100[i,which(gmpg100[i,]<=thres)])
}

qpg100 <- graphNEL(nodes=my.nodes.pg100, edgeL=edLpg100)

tiff("figure9_4.tif",width=5,height=5,units="in",res=1200)
par(mfrow=c(2,1))
qpPlotNetwork(qpg100,minimumSizeConnComp=4)
mtext("A", side = 3, line = 0, adj = 0.01, cex = 1, padj = 1)
qpPlotNetwork(pcc100,minimumSizeConnComp=4)
mtext("B", side = 3, line = 0, adj = 0.01, cex = 1)
#par(mfrow=c(1,1))
dev.off()

tiff("figure9_5.tif",width=5,height=5,units="in",res=1200)
par(mfrow=c(2,1))
qpPlotNetwork(qpg100,vertexSubset="Lhx2", boundary=TRUE)
mtext("A", side = 3, line = 0, adj = 0.01, cex = 1, padj = 1)
qpPlotNetwork(pcc100,vertexSubset="Lhx2", boundary=TRUE)
mtext("B", side = 3, line = 0, adj = 0.01, cex = 1)
dev.off()

rm(list = ls())

#################################################################
# Master regulators

# figure 9.6

library(Biobase)
load("eset_tf.rda")

tiff("figure9_6.tif",width=5,height=5,units="in",res=1200)
par(mfrow=c(1,2),mai=c(0.52,0.52,0.72,0.12))
plot(density(exprs(eset.tf)),
     main="Raw gene expression",xlab="",ylab="",cex.main=0.9)
mtext("A", side = 3, line = 0, adj = 0, cex = 1)
plot(density(log2(exprs(eset.tf)),na.rm=T),
     main="Logarithm\nof gene expression",xlab="",ylab="",cex.main=0.9)
mtext("B", side = 3, line = 0, adj = 0, cex = 1)
dev.off()

# figure 9.7
library(RTN)
annot2<-read.csv("annot.csv",row.names=1)

target.tf.probes<-c("1341_at","40511_at","41504_s_at","33592_at")
names(target.tf.probes)<-c("SPI1","GATA3","MAF","ZBTB7B")

tf.rtni<-new("TNI", 
             gexp=exprs(eset.tf), 
             transcriptionFactors=target.tf.probes)
tf.rtni<-tni.preprocess(tf.rtni, gexpIDs=annot2)
tf.rtni<-tni.permutation(tf.rtni, estimator='kendall', pValueCutoff=0.03)
tf.rtni<-tni.bootstrap(tf.rtni, estimator='kendall', consensus=95)
tf.rtni<-tni.dpi.filter(tf.rtni)

tiff("figure9_7.tif",width=5,height=5,units="in",res=1200)
par(mfrow=c(1,2),mai=c(0.52,0.52,0.72,0.12))
g<-tni.graph(tf.rtni)
#V(g)$color <- "black"
plot(g,vertex.shape="none",vertex.label.color="black",vertex.label.cex=0.5)
mtext("A", side = 3, line = 0, adj = 0.5, cex = 1)
V(g)$label<-as.character(annot2[annot2$PROBEID %in% V(g)$name, "SYMBOL"])
plot(g,vertex.shape="none",vertex.label.color="black",vertex.label.cex=0.5)
mtext("B", side = 3, line = 0, adj = 0.5, cex = 1)
dev.off()

rm(list = ls())

#################################################################
# GeneAnswers demo

library(GeneAnswers)
furin.genes<-read.csv("furin_significant_gids.csv")
names(furin.genes)[1]<-"GeneID"
furin.input<-data.frame("Entrez Gene ID"=furin.genes$Gene, fold.change=log2(furin.genes$Activated.KO/furin.genes$Activated.WT))
genAns<-geneAnswersBuilder(furin.input, 
                           'org.Mm.eg.db', 
                           categoryType='KEGG', 
                           known=T, 
                           geneExpressionProfile=furin.genes)
genAnsRead<-geneAnswersReadable(genAns)

geneAnswersChartPlots(genAnsRead, 
                      chartType='pieChart',
                      newWindow=F,cex=0.6)

geneAnswersConceptNet(genAnsRead, 
                      colorValueColumn='fold.change', 
                      centroidSize='pvalue',output='interactive')

geneAnswersConceptNet(genAnsRead,
                      colorValueColumn='fold.change',
                      centroidSize='geneNum',output='interactive')

tiff("figure9_9_color.tif",width=5,height=5,units="in",res=1200)
geneAnswersHeatmap(genAns, catTerm=TRUE, geneSymbol=TRUE)
dev.off()
