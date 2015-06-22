library(topGO)
data(geneList)

library(hgu95av2.db)

affyLib <-"hgu95av2.db"

#smt<-function(p){return(p<0.001)}

GOdata <- new("topGOdata", ontology = "BP",allGenes = geneList, geneSel = function(p){return(p<0.001)}, nodeSize = 10, annot = annFUN.db, affyLib = "hgu95av2.db")

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "elim", statistic = "KS")

allRes.BP <- GenTable(GOdata, classicKS = resultFisher, elimKS = resultKS, orderBy = "elimKS", topNodes=20)

showSigOfNodes(GOdata, score(resultKS), firstSigNodes = 5, useInfo = 'all')

