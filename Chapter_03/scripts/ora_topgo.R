immunome<-read.csv("human_gene_set.txt")
target_set<-read.csv("human_gene_selected_set.txt")


library(org.Hs.eg.db)
goannot2<-select(org.Hs.eg.db,keys=immunome$EntrezGeneID,cols="GO",keytype="ENTREZID")

# Prepare the gene GO annotation in proper format into gene2GO

a <- 1
gene2GO<-c()
genevec <- c()

for (genid in immunome$EntrezGeneID){

  genevec[a] <- genid
  gene2GO[a] <- list(goannot2[goannot2$ENTREZID==genid,"GO"])
  a <- a+1
}
names(gene2GO)<-genevec

# create geneList to mark the target genes in a 0,1 vector

geneList<-rep(0,length(immunome$EntrezGeneID))
geneList[immunome$EntrezGeneID %in% target_set$EntrezGeneID]<-1
names(geneList)<-genevec
geneList<-factor(geneList)

library(topGO)

GOdata.MF <- new("topGOdata", ontology = "MF", description="MF on innate immunity genes",allGenes = geneList, annot =annFUN.gene2GO, gene2GO = gene2GO)

resultFisher.MF.classic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
allRes.MF.classic <- GenTable(GOdata, classicFisher = resultFisher.MF.classic,topNodes=20)

resultFisher.MF.weight <- runTest(GOdata, algorithm = "weight", statistic = "fisher")
allRes.MF.weight <- GenTable(GOdata, classicFisher = resultFisher.MF.weight,topNodes=20)

resultFisher.MF.parentchild <- runTest(GOdata, algorithm = "parentchild", statistic = "fisher")
allRes.MF.parentchild <- GenTable(GOdata, classicFisher = resultFisher.MF.parentchild,topNodes=20)

