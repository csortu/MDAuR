immunome<-read.csv("human_gene_set.txt")
target_set<-read.csv("human_gene_selected_set.txt")


library("GOstats")
library("org.Hs.eg.db")

# GO analysis
paramsGO <- new("GOHyperGParams", geneIds = target_set$EntrezGeneID, universeGeneIds =immunome$EntrezGeneID, annotation = "org.Hs.eg.db", ontology = "BP", pvalueCutoff = 0.01, conditional = FALSE, testDirection = "over")

Over.GO <- hyperGTest(paramsGO)
head(summary(Over.GO))

# KEGG analysis

paramsKEGG <- new("KEGGHyperGParams", geneIds = target_set$EntrezGeneID, universeGeneIds =immunome$EntrezGeneID, annotation = "org.Hs.eg.db", pvalueCutoff = 0.01, testDirection = "over")

Over.KEGG <- hyperGTest(paramsKEGG)
head(summary(Over.KEGG))



