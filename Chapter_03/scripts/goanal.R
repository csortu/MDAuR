library(topGO)
vulfile <- "/home/ortutay/Dokumentumok/Sync/ImmuneGenes/vul.csv"
gene2gofile <- "/home/ortutay/Dokumentumok/Sync/ImmuneGenes/gene2go.csv"

vul <- read.table(vulfile)
genecats<-read.table("/home/ortutay/Dokumentumok/Sync/ImmuneGenes/genecats.csv")

lines<-readLines(gene2gofile)

a <- 1
gene2GO<-c()
genevec <- c()

for (row in lines){

  rowvec <- unlist(strsplit(row," "))
  genevec[a] <- rowvec[1]
  gene2GO[a] <- list(rowvec[2:length(rowvec)])

  a <- a+1
}
names(gene2GO)<-genevec

genList<-rep(0,length(genevec))
names(genList)<-genevec

for (a in 1:length(levels(genecats[,1]))){
  type<-levels(genecats[,1])[a]

  for (gene in as(genecats[genecats[,1]==type,3],"character")){
    if(is.na(genList[gene])){
    }else{
      genList[gene]<-1
    }
  }
}

GOdata <- new("topGOdata", ontology = "MF", allGenes = genList, annot =annFUN.gene2GO, gene2GO = gene2GO)

test.classic.stat <- new("classicCount", testStatistic = GOFisherTest,name = "Fisher test")
resultFis <- getSigGroups(GOdata, test.classic.stat)
printGraph(GOdata, resultFis, firstSigNodes = 9, fn.prefix = "tGO",pdfSW = TRUE)

test.elim.stat <- new("elimCount", testStatistic = GOFisherTest,name = "Fisher test", cutOff = 0.01)
resultElim <- getSigGroups(GOdata, test.elim.stat)
printGraph(GOdata, resultElim, firstSigNodes = 9, fn.prefix = "tGO",pdfSW = TRUE)

test.weight.stat <- new("weightCount", testStatistic = GOFisherTest,name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.weight.stat)
printGraph(GOdata, resultWeight, firstSigNodes = 9, fn.prefix = "tGO",pdfSW = TRUE)

l <- list(classic = score(resultFis),elim = score(resultElim), weight = score(resultWeight))
allRes <- genTable(GOdata, l, orderBy = "weight", ranksOf = "classic",top = 20)

# The vulnerability score test
# Testing the 50 genes with the highest vulnerability

genList<-rep(0,length(genevec))
names(genList)<-genevec

for (gene in as(vul$V1[vul$V2>0.005],"character")){
  genList[gene]<-1
}


