library(topGO)
grstatfile <- "/home/ortutay/Dokumentumok/Sync/ImmuneGenes/graphanal.csv"
gene2gofile <- "/home/ortutay/Dokumentumok/Sync/ImmuneGenes/gene2go.csv"
genecatfile <- "/home/ortutay/Dokumentumok/Sync/ImmuneGenes/genecats.csv"

genecats<-read.table(genecatfile)
grstat <- read.table(grstatfile)
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

# GO vs disease categories analysis

# Create genList which is a named list of 0s and 1s
# 1 shows the selected genes

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

genList<-factor(genList)


for (gotype in c("MF","BP","CC")){

  GOdata <- new("topGOdata", ontology = gotype, allGenes = genList, annot =annFUN.gene2GO, gene2GO = gene2GO)

# Calssic test
  
  test.classic.stat <- new("classicCount", testStatistic = GOFisherTest,name = "Fisher test")
  resultFis <- getSigGroups(GOdata, test.classic.stat)
  printGraph(GOdata, resultFis, firstSigNodes = 9, fn.prefix = paste(gotype,"ontology_disease",sep=""),pdfSW = TRUE)

# elim test

  test.elim.stat <- new("elimCount", testStatistic = GOFisherTest,name = "Fisher test", cutOff = 0.01)
  resultElim <- getSigGroups(GOdata, test.elim.stat)
  printGraph(GOdata, resultElim, firstSigNodes = 9, fn.prefix = paste(gotype,"ontology_disease",sep=""),pdfSW = TRUE)

# weight test

  test.weight.stat <- new("weightCount", testStatistic = GOFisherTest,name = "Fisher test", sigRatio = "ratio")
  resultWeight <- getSigGroups(GOdata, test.weight.stat)
  printGraph(GOdata, resultWeight, firstSigNodes = 9, fn.prefix = paste(gotype,"ontology_disease",sep=""),pdfSW = TRUE)

  # Make the table from the results
  
  l <- list(classic = score(resultFis),elim = score(resultElim), weight = score(resultWeight))
  allRes <- genTable(GOdata, l, orderBy = "weight", ranksOf = "classic",top = 20)
  
  write.table(allRes,file=paste(gotype,"ontology_disease_results.csv",sep=""),sep = "\t") 
}

# The score tests
# Testing the 50 genes with the highest scores


# Find out the top 50 genes limits for the categories

limit<-c()
limit["deg"]<-sort(grstat$deg,decreasing=T)[50]
limit["vul"]<-sort(grstat$vul,decreasing=T)[50]
limit["clos"]<-sort(grstat$clos,decreasing=T)[50]

for (sctype in c("deg","vul","clos")){
  genList<-rep(0,length(genevec))
  names(genList)<-genevec
  
  for (gene in row.names(grstat[grstat[sctype]>limit[sctype],])){
    if(is.na(genList[gene])){
    }else{
      genList[gene]<-1
    }
  }
  genList<-factor(genList)
  for (gotype in c("MF","BP","CC")){

    GOdata <- new("topGOdata", ontology = gotype, allGenes = genList, annot =annFUN.gene2GO, gene2GO = gene2GO)
    
# Calssic test
  
    test.classic.stat <- new("classicCount", testStatistic = GOFisherTest,name = "Fisher test")
    resultFis <- getSigGroups(GOdata, test.classic.stat)
    printGraph(GOdata, resultFis, firstSigNodes = 9, fn.prefix = paste(gotype,sctype,"_ontology_score",sep="_"),pdfSW = TRUE)

# elim test

    test.elim.stat <- new("elimCount", testStatistic = GOFisherTest,name = "Fisher test", cutOff = 0.01)
    resultElim <- getSigGroups(GOdata, test.elim.stat)
    printGraph(GOdata, resultElim, firstSigNodes = 9, fn.prefix = paste(gotype,sctype,"_ontology_score",sep="_"),pdfSW = TRUE)

# weight test

    test.weight.stat <- new("weightCount", testStatistic = GOFisherTest,name = "Fisher test", sigRatio = "ratio")
    resultWeight <- getSigGroups(GOdata, test.weight.stat)
    printGraph(GOdata, resultWeight, firstSigNodes = 9, fn.prefix = paste(gotype,sctype,"_ontology_score",sep="_"),pdfSW = TRUE)

  # Make the table from the results
  
    l <- list(classic = score(resultFis),elim = score(resultElim), weight = score(resultWeight))
    allRes <- genTable(GOdata, l, orderBy = "weight", ranksOf = "classic",top = 20)
    
    write.table(allRes,file=paste(gotype,sctype,"ontology_score_results.csv",sep="_"),sep = "\t") 
  }

}

