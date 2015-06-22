mouse_genes<-read.csv("mouse_ranked_gene_set.txt")

library(org.Mm.eg.db)

pathannot<-select(org.Mm.eg.db,keys=mouse_genes$EntrezGeneID,cols="PATH",keytype="ENTREZID")
pathannot<-pathannot[complete.cases(pathannot),]

bigpaths<-names(table(pathannot$PATH)[table(pathannot$PATH)>100])


library(KEGG.db)
keggid2keggname <- as.list(KEGGPATHID2NAME)
#keggid2keggname[bigpaths]

for (pathid in bigpaths){
  pathgenes<-pathannot$ENTREZID[pathannot$PATH==pathid]
  gene.set<-as.numeric(mouse_genes$EntrezGeneID %in% pathgenes)
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

}
	
for (pathid in bigpaths){
  pathgenes<-pathannot$ENTREZID[pathannot$PATH==pathid]
  gene.set<-as.numeric(mouse_genes$EntrezGeneID %in% pathgenes)
  not.gene.set<-1-gene.set
  scale.set    <- sqrt((length(gene.set) - sum(gene.set))/sum(gene.set))
  scale.not.set <- sqrt(sum(gene.set)/(length(gene.set) - sum(gene.set)))
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

# permutation test

  groupsize<-sum(gene.set)
  esdist<-c()

  for (i in 1:1000){
    rnd.gene.set<-as.numeric(mouse_genes$EntrezGeneID %in% sample(mouse_genes$EntrezGeneID,groupsize))
    not.gene.set<-1-rnd.gene.set
    RRES<-cumsum(gene.set*scale.set-not.gene.set*scale.not.set)
    max.ES <- max(RRES)
    min.ES <- min(RRES)
    if (max.ES > - min.ES) {
       ES.random <- signif(max.ES, digits=5)
    } else {
       ES.random <- signif(min.ES, digits=5)
    }

    esdist[i]<-abs(ES.random)
  }

  ttres<-t.test(esdist,alternative="less",mu=abs(ES))

  jfile<-paste(getwd(),.Platform$file.sep,pathid,'.jpg',sep="")
  jpeg(filename=jfile)
  plot(RES,type="l",main=unlist(keggid2keggname[pathid]),sub=paste("p =",ttres$p.value))
  text(arg.ES,ES,label=ES)
  text(arg.ES,ttres$estimate,label=paste("Random ES: ",ttres$estimate))
  dev.off()
}


OLD.GSEA.EnrichmentScore <- function(gene.list, gene.set) {  
#
# Computes the original GSEA score from Mootha et al 2003 of gene.set in gene.list 
#
# Inputs:
#   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
#   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
#
# Outputs:
#   ES: Enrichment score (real number between -1 and +1) 
#   arg.ES: Location in gene.list where the peak running enrichment occurs (peak of the "mountain") 
#   RES: Numerical vector containing the running enrichment score for all locations in the gene list 
#   tag.indicator: Binary vector indicating the location of the gene sets (1's) in the gene list 
#

   tag.indicator <- sign(match(gene.list, gene.set, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
   no.tag.indicator <- 1 - tag.indicator 
   N <- length(gene.list) 
   Nh <- length(gene.set) 
   Nm <-  N - Nh 

   norm.tag    <- sqrt((N - Nh)/Nh)
   norm.no.tag <- sqrt(Nh/(N - Nh))

   RES <- cumsum(tag.indicator * norm.tag - no.tag.indicator * norm.no.tag)      
   max.ES <- max(RES)
   min.ES <- min(RES)
   if (max.ES > - min.ES) {
      ES <- signif(max.ES, digits=5)
      arg.ES <- which.max(RES)
   } else {
      ES <- signif(min.ES, digits=5)
      arg.ES <- which.min(RES)
   }
   return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))    
}

# permutation test
groupsize<-sum(gene.set)
scale.set    <- sqrt((length(gene.set) - sum(gene.set))/sum(gene.set))
scale.not.set <- sqrt(sum(gene.set)/(length(gene.set) - sum(gene.set)))
esdist<-c()

for (i in 1:1000){
  rnd.gene.set<-as.numeric(mouse_genes$EntrezGeneID %in% sample(mouse_genes$EntrezGeneID,groupsize))
  not.gene.set<-1-rnd.gene.set
  RRES<-cumsum(gene.set*scale.set-not.gene.set*scale.not.set)
  max.ES <- max(RRES)
  min.ES <- min(RRES)
  if (max.ES > - min.ES) {
     ES <- signif(max.ES, digits=5)
  } else {
     ES <- signif(min.ES, digits=5)
  }

  esdist[i]<-abs(ES)
}
plot(density(esdist))
