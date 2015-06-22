library(Cairo)

col <- gray.colors(4)

library(seqinr)
ttn<-read.fasta(file="human_ttn_genomic_dna.fasta")
seq<-ttn[[1]]
seq.rho<-rho(seq) # the rho statistics
seq.zsq.base<-zscore(seq, model = "base", exact = F) # using the 'base' model
seq.zsq.codon<-zscore(seq, model = "codon", exact = F) # using the 'codon' model

opar <- par()

## CairoTIFF("figure2_1.tif",width=5,height=5,units="in",res=1200)

## lwd <- 4

## par(mfrow=c(2,2),mai=c(0.42,0.67,0.42,0.12))
## plot(seq.rho,ylim=c(min(seq.rho),max(seq.rho)), lwd = lwd, col=col, main="Rho")
## mtext("A", side = 3, line = -0.2, adj = -0.3, cex = 1.3)


## plot(seq.rho-1,ylim=c(min(seq.rho)-1,max(seq.rho)-1), lwd = lwd, col=col, main="Rho – 1")
## mtext("B", side = 3, line = -0.2, adj = -0.3, cex = 1.3)

## plot(seq.zsq.codon,ylim=c(min(seq.zsq.codon),max(seq.zsq.codon)), lwd = lwd, col=col, main="Z-score\nCodon model")
## mtext("C", side = 3, line = -0.2, adj = -0.3, cex = 1.3)

## plot(seq.zsq.base,ylim=c(min(seq.zsq.base),max(seq.zsq.base)), lwd = lwd, col=col, main="Z-score\nBase model")
## mtext("D", side = 3, line = -0.2, adj = -0.3, cex = 1.3)

## dev.off()

CairoTIFF("figure2_2.tif",width=5,height=5,units="in",res=1200)

par(opar)

par(mfrow=c(2,2),mai=c(0.22,0.67,0.42,0.12),xaxt="n",yaxt="n",mgp=c(2,1,0))
for (f.i in 1:4){
  win.size <- c(1000,5000,10000,50000)[f.i]
  gc.perc<-vector()
  k<-1
  for (i in seq(from=1,to=length(seq)-win.size,by=win.size/10)){
    j<-i+win.size-1
    gc.perc[k]<-GC(seq[i:j])
    k<-k+1
  }
  plot(gc.perc*100,type="l", main=paste("Window size",win.size),ylab="GC%",xlab="",col.axis="white")
  mtext(LETTERS[f.i], side = 3, line = -0.2, adj = -0.3, cex = 1.3)
  axis(1,labels=F)
  axis(2,labels=F)
}

dev.off()

## rhod<-read.fasta("human_rhodopsin_protein.fasta",seqtype = "AA")
## rhodseq<-rhod[[1]]

## CairoTIFF("figure2_3.tif",width=5,height=5,units="in",res=1200)
## AAstat(rhodseq)
## dev.off()

## CairoTIFF("figure2_4.tif",width=5,height=5,units="in",res=1200)
## hpat<-read.table("aa_hydropathy.txt")
## hpatseq<-hpat[rhodseq[1:length(rhodseq)],2]
## par(mfrow=c(1,1),mai=c(0.82,0.82,0.42,0.12))
## for (f.i in 4:4){
##   win.size <- c(2,4,8,16)[f.i]
##   hpat.m<-vector()
##   k<-1
##   for (i in seq(from=1,to=length(hpatseq)-win.size)){
##     j<-i+win.size-1
##     hpat.m[k]<-mean(hpatseq[i:j])
##     k<-k+1
##   }
##   plot(hpat.m,type="l", main="Rhodopsin protein", ylab="Hydropathy", xlab="Sequence position")
##   mtext(LETTERS[f.i], side = 3, line = -0.2, adj = -0.3, cex = 1.3)
##   arrows(x0=c(10,38,75,125,155,205,275,288),x1=c(10,38,75,125,155,205,258,288),y0=c(1.2,2.7,2.7,2.7,2.7,2.7,2.7,2),y1=c(0.65,2.4,2.1,2,2.2,2.4,2.6,1.5),code=2,length=0.1)
## }

## dev.off()


## library(GenomeGraphs)
## mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

## gene <- makeGene(id = "ENSG00000140564",type = "ensembl_gene_id", biomart = mart, dp=DisplayPars(color="black"))

## transcript <- makeTranscript(id = "ENSG00000140564",type = "ensembl_gene_id", biomart = mart, dp=DisplayPars(color="gray"))

## customann<-makeAnnotationTrack(start=c(90876557,90876757,90878557), end=c(90876597,90876897,90878697), feature=c('bind','bind','del'), dp=DisplayPars(bind = 'gray', del='black'))

## my.gene<-getBM(c('chromosome_name', 'start_position', 'end_position'), filters=c('ensembl_gene_id'), value="ENSG00000140564", mart = mart)

## snpmart <- useMart("snp",dataset="hsapiens_snp")
## snps<-getBM(c('refsnp_id', 'allele', 'chrom_start'), filters = c('chr_name', 'start', 'end'), values= list(my.gene$chromosome_name, my.gene$start_position, my.gene$end_position), mart = snpmart)

## snpannot<-makeAnnotationTrack(start=snps$chrom_start, end=snps$chrom_start)


## CairoTIFF("figure2_5.tif",width=5,height=5,units="in",res=1200)

## gdPlot(list(makeTitle("Human furin gene"), makeIdeogram(chromosome = 15), "FURIN" = gene, "RNA transcripts" = transcript, "SNP" = snpannot,"Reg." = customann,makeGenomeAxis(dp = DisplayPars(cex = 0.6))))

## dev.off()

