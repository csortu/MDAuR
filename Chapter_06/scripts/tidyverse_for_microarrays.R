library(GEOquery)
library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)

###########################################################################
#
# This script demonstrates how to use tidyverse functions to easily 
# manipulate microarray data
#
###########################################################################


## Some basic tidying work on example dataframe

example.df <- data.frame(Sample1=rnorm(5,mean = 4,sd = 1),
                         Sample2=rnorm(5,mean = 4,sd = 1),
                         Sample3=rnorm(5,mean = 4,sd = 1),
                         row.names = paste0("probe",1:5))

example.df.mutated <- example.df %>% mutate(ID=rownames(example.df))
tidy.df <- example.df.mutated %>% gather("Chip","Expr",1:3)

my.plot <- ggplot(data=tidy.df,aes(x=Chip,y=Expr))
my.plot + geom_point()
my.plot + geom_line()
my.plot + geom_boxplot()

## Let's do the same on real gene expression data

datadir <- getwd() # You will need to work in a directory with at least 80 MB free space


gse <- getGEO("GSE43805",destdir=datadir)
gse
gse[[1]]
microarray.dat <- as.data.frame(exprs(gse[[1]]))

microarray.tidy <- microarray.dat %>% 
  mutate(ID=rownames(microarray.dat)) %>%
  gather("Chip","Expr",1:4)

ggplot(data=microarray.tidy,aes(x=Chip,y=Expr)) + geom_boxplot()
ggplot(data=microarray.tidy,aes(x=Chip,y=log2(Expr))) + geom_boxplot()

## Add annotation with join

example.annot <- data.frame(probeID=paste0("probe",1:5),
                            Symbol=c(rep("CCT3",3),
                                     "IL6","BTK"),
                            GeneID = c(rep(7203,3),
                                       3569,695))
left_join(tidy.df,example.annot,by=c("ID"="probeID"))
tidy.df %>% left_join(example.annot,by=c("ID"="probeID"))
tidy.df %<>% left_join(example.annot,by=c("ID"="probeID"))

# Real annottion

gpl <- getGEO("GPL10558",destdir=datadir)
gpl

Meta(gpl)$title
colnames(Table(gpl))
Table(gpl)[1001:1005,1:18]

annot.dat <- Table(gpl)

microarray.tidy %<>% left_join(annot.dat,by=c("ID"))

# Data preparation is done, let's do basic calculations and visualization

gene.expr <- tidy.df %>% group_by(Chip,GeneID,Symbol) %>%
  summarise(MeanExpr=mean(Expr))

gene.expr

ggplot(data=gene.expr,aes(x=Chip,y=MeanExpr,color=Symbol,group=GeneID)) + 
  geom_point() + 
  geom_line()

# Now the same for real data

my.genes <- c("CCT3","IL6","RGL1","HBE1")

selected.genes <- microarray.tidy %>% group_by(Chip,Entrez_Gene_ID,Symbol) %>%
  filter(Symbol %in% my.genes) %>%
  summarise(MeanExpr=mean(Expr))

ggplot(data=selected.genes,
       aes(x=Chip,y=MeanExpr,color=Symbol,group=Entrez_Gene_ID)) + 
  geom_point() + geom_line() + 
  ggtitle("Expression of selected genes in TMEM88 KO stem cells") + 
  annotate("text",x=c(1.5,3.5),y=c(8.5,10),
           label=c("Control\nsamples","TMEM88 KO")) + 
  annotate("segment",x=2.5,xend = 2.5,
           y=6.9,yend=11.3,color="gray")

