library( GenomeGraphs)

#listMarts()

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#listDatasets(mart)

gene <- makeGene(id = "ENSG00000140564",type = "ensembl_gene_id", biomart = mart)
gdPlot(gene)

transcript <- makeTranscript(id = "ENSG00000140564",type = "ensembl_gene_id", biomart = mart)
gdPlot(list(gene, transcript))

gdPlot(list(makeTitle("Human furin gene"), makeIdeogram(chromosome = 15), gene, transcript, makeGenomeAxis()))


# get gene start and end pos

my.gene<-getBM(c('chromosome_name','start_position','end_position'),filters = c('ensembl_gene_id'),value = "ENSG00000140564", mart = mart)

# get snp tranck on the plot

snpmart <- useMart("snp",dataset="hsapiens_snp")


#snps<-getBM(c('refsnp_id','allele','chrom_start','clinical_significance'),filters = c('chr_name','chrom_start','chrom_end'), values = list(my.gene$chromosome_name,my.gene$start_position,my.gene$end_position), mart = snpmart)

snps<-getBM(c('refsnp_id','allele','chrom_start'),filters = c('chr_name','start','end'), values = list(my.gene$chromosome_name,my.gene$start_position,my.gene$end_position), mart = snpmart)


snpannot<-makeAnnotationTrack(start=snps$chrom_start,end=snps$chrom_start)


gdPlot(list(makeTitle("Human furin gene"), makeIdeogram(chromosome = 15), gene, transcript, snpannot,makeGenomeAxis()))

# create custom annotation

customann<-makeAnnotationTrack(start=c(90876557,90876757,90878557),end=c(90876597,90876897,90878697),feature=c('bind','bind','del'),dp = DisplayPars(bind = 'blue',del='red'))


gdPlot(list(makeTitle("Human furin gene"), makeIdeogram(chromosome = 15), gene, transcript, snpannot,customann,makeGenomeAxis()))
