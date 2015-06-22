library(biomaRt)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

immunome<-read.csv("human_gene_set.txt")

#goannot<-getBM(c("entrezgene","go_id"),filters="entrezgene",values=c(963,965),mart=mart)

goannot<-getBM(c("entrezgene","go_id","name_1006"),filters="entrezgene",values=immunome,mart=mart)

goannot[goannot$entrezgene==965,]

# Do the same with org.Hs.eg.db

library(org.Hs.eg.db)

immunome<-read.csv("human_gene_set.txt")

goannot2<-select(org.Hs.eg.db,keys=as.character(immunome$EntrezGeneID),columns="GO",keytype="ENTREZID")

pathannot<-select(org.Hs.eg.db,keys=as.character(immunome$EntrezGeneID),columns="PATH",keytype="ENTREZID")

clannot<-select(org.Hs.eg.db,keys=as.character(immunome$EntrezGeneID),columns="CHRLOC",keytype="ENTREZID")

keytypes(org.Hs.eg.db)
