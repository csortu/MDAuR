my.ala <- read.table("mouse_ala.txt")$V1
my.tca <- read.table("mouse_tca.txt")$V1
my.gly <- read.table("mouse_glyco.txt")$V1


library(org.Mm.eg.db)

mouse.keys<-keys(org.Mm.eg.db,keytype="ENTREZID")
pathannot<-select(org.Mm.eg.db,keys=mouse.keys,columns="PATH",keytype="ENTREZID")
pathannot<-pathannot[complete.cases(pathannot),]
my.path <- "00250"


get.RES <- function(mouse.genes,pathgenes){

  gene.set<-as.numeric(mouse.genes %in% pathgenes)
  not.gene.set<-1-gene.set
  scale<-sum(gene.set)/sum(not.gene.set)
  scale.set    <- sqrt((length(gene.set) - sum(gene.set))/sum(gene.set))
  scale.not.set <- sqrt(sum(gene.set)/(length(gene.set) - sum(gene.set)))
  RES.simple<-cumsum(gene.set-not.gene.set*scale)
  RES<-cumsum(gene.set*scale.set-not.gene.set*scale.not.set)
  max.ES <- max(RES)
  min.ES <- min(RES)
  if (max.ES > - min.ES) {
    ES <- signif(max.ES, digits=5)
    arg.ES <- which.max(RES)
  } else {
    ES <- signif(min.ES, digits=5)
    arg.ES <- which.min(RES)
  }

  return(list(RES,ES))
}
	
pathgenes<-pathannot$ENTREZID[pathannot$PATH==my.path]

library(Cairo)



CairoTIFF("figure3_2.tif",width=5,height=5,units="in",res=1200)

#par(mfrow=c(2,2),mai=c(0.22,0.67,0.42,0.12),xaxt="n",yaxt="n",mgp=c(2,1,0))
par(mfrow=c(2,2),mai=c(0.42,0.52,0.42,0.12))

mouse.genes <- c(sample(c(my.ala,my.tca)),sample(c(my.gly,my.gly,my.gly,my.gly)))
plot(get.RES(mouse.genes,pathgenes)[[1]],type="l",xlab="",ylab="",main="Top genes")
mtext("A", side = 3, line = 0.6, adj = -0.25, cex = 1)


mouse.genes <- c(sample(c(my.gly,my.gly,my.gly,my.gly)),sample(c(my.ala,my.tca)))
plot(get.RES(mouse.genes,pathgenes)[[1]],type="l",xlab="",ylab="",main="Bottom genes")
mtext("B", side = 3, line = 0.6, adj = -0.25, cex = 1)


mouse.genes <- c(sample(c(my.gly,my.gly)),sample(c(my.ala,my.tca)),sample(c(my.gly,my.gly)))
plot(get.RES(mouse.genes,pathgenes)[[1]],type="l",xlab="",ylab="",main="Middle genes")
mtext("C", side = 3, line = 0.6, adj = -0.25, cex = 1)
abline(h=0,lwd=0.5)

mouse.genes <- sample(c(my.gly,my.gly,my.gly,my.ala,my.tca))
plot(get.RES(mouse.genes,pathgenes)[[1]],type="l",xlab="",ylab="",main="No association")
mtext("D", side = 3, line = 0.6, adj = -0.25, cex = 1)
abline(h=0,lwd=0.5)



dev.off()

mouse.genes<-read.csv("mouse_ranked_gene_set.txt")
pathannot<-select(org.Mm.eg.db, keys=as.character(mouse_genes$EntrezGeneID), columns="PATH",keytype="ENTREZID")
pathannot<-pathannot[complete.cases(pathannot),]

pathid<-"05200"
pathgenes<-pathannot$ENTREZID[pathannot$PATH==pathid]
gene.set<-as.numeric(mouse_genes$EntrezGeneID %in% pathgenes)
gene.set-not.gene.set

scale.set <- sqrt((length(gene.set) - sum(gene.set))/sum(gene.set))
scale.not.set <- sqrt(sum(gene.set)/(length(gene.set) - sum(gene.set)))

RES<-cumsum(gene.set*scale.set-not.gene.set*scale.not.set)

CairoTIFF("figure3_2.tif",width=5,height=5,units="in",res=1200)

plot(RES,type="l",main="Pathways in cancer",xlab="Genes in order")

dev.off()

