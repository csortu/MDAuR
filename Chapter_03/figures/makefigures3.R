library(topGO)
data(geneList)
library(hgu95av2.db)
library(Cairo)


GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, geneSel = function(p){return(p<0.01)}, nodeSize = 10, annot = annFUN.db, affyLib = "hgu95av2.db")

##resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")

resultKS <- runTest(GOdata, algorithm = "elim", statistic = "KS")

CairoTIFF("figure3_1.tif",width=5,height=5,units="in",res=1200)

showSigOfNodes(GOdata, score(resultKS), firstSigNodes = 2)

dev.off()

