library(GeneAnswers)

furin.genes<-read.csv("furin_significant_gids.csv")
names(furin.genes)[1]<-"GeneID"

furin.input<-data.frame("Entrez Gene ID"=furin.genes$Gene,fold.change=log2(furin.genes$Activated.KO/furin.genes$Activated.WT))

genAns<-geneAnswersBuilder(furin.input, 'org.Mm.eg.db', categoryType='KEGG', known=T, geneExpressionProfile=furin.genes)
genAnsRead<-geneAnswersReadable(genAns)

geneAnswersChartPlots(genAnsRead, chartType='all')

geneAnswersConceptNet(genAnsRead,colorValueColumn='fold.change', centroidSize='pvalue',output='fixed')
geneAnswersConceptNet(genAnsRead,colorValueColumn='fold.change', centroidSize='geneNum',output='fixed')
geneAnswersConceptNet(genAnsRead,colorValueColumn='fold.change', centroidSize='pvalue',output='interactive')

geneAnswersHeatmap(genAns, catTerm=TRUE, geneSymbol=TRUE)

topPATH(genAns,keepID=T)
#topPATH(genAns)



#genAnsBP<-geneAnswersBuilder(furin.input, 'org.Mm.eg.db', categoryType='GO.BP', known=T, FDR.correction=TRUE,pvalueT=0.01)
#genAnsReadBP<-geneAnswersReadable(genAnsBP)

#geneAnswersConceptNet(genAnsReadBP,colorValueColumn='fold.change', centroidSize='pvalue',output='fixed')
#geneAnswersConceptNet(genAnsReadBP,colorValueColumn='fold.change', centroidSize='pvalue',output='interactive')



#genAnsMF<-geneAnswersBuilder(furin.input, 'org.Mm.eg.db', categoryType='GO.MF', known=T, FDR.correction=TRUE,pvalueT=0.01)
#genAnsReadMF<-geneAnswersReadable(genAnsMF)

# interest genes:
# positives:
# Fgf2, Lama1, Mapk10, Pik3cd, Egfr, Vegfa, H2-Ob, H2-b1
# IDs: 14173, 16772, 26414, 18707, 13649, 22339, 15002, 14963
# negatives:
# Ifng, Mapk8, Il3, Ccl3, Csf2, Cd40, Tlr4, Cxcl13
# IDs: 15978, 26419, 16187, 20302, 12981, 21939, 21898, 55985

#target.genes.up<-c(14173, 16772, 26414, 18707, 13649, 22339, 15002, 14963)
#target.genes.down<-c(15978, 26419, 16187, 20302, 12981, 21939, 21898, 55985)

