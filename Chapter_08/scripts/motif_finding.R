bed.file <- 'GSM1504361_PU.1_ChIPseq_PU.1_B_cells.bed.gz'
peak.tab <- read.delim(bed.file,header=F)
names(peak.tab) <- c('chr','start','end','peakName','score')
dim(peak.tab)

library(rGADEM)

peak.ranges <- RangedData(IRanges(start=peak.tab$start,end=peak.tab$end),space=paste0('chr',peak.tab$chr))


# annotation!
# MGSCv37 is the same as UCSC mm9! find correct bioconductor package
# https://genome.ucsc.edu/FAQ/FAQreleases.html


library(BSgenome.Mmusculus.UCSC.mm9)

gadem <- GADEM(peak.ranges,verbose=1,genome=Mmusculus) # takes long-long time 3 hours on my computer

gadem

nMotifs(gadem)

nOccurrences(gadem)

consensus(gadem)

bestmotif <- which(nOccurrences(gadem)==max(nOccurrences(gadem)))

consensus(gadem)[bestmotif]

getPWM(gadem)
motifs <- getPWM(gadem)

seqLogo(motifs[[bestmotif]])
seqLogo(reverseComplement(motifs[[bestmotif]]))


save(gadem,file="pu1_gadem.rda")
load("pu1_gadem.rda")


############ run only on a subset of peaks
## 50 best

peak.tab.sel <- peak.tab[peak.tab$score %in% tail(sort(peak.tab$score),50),]

peak.ranges.sel <- RangedData(IRanges(start=peak.tab.sel$start,end=peak.tab.sel$end),space=paste0('chr',peak.tab.sel$chr))

gadem.sel <- GADEM(peak.ranges.sel,verbose=1,genome=Mmusculus)

gadem.sel

nMotifs(gadem.sel)

nOccurrences(gadem.sel)

consensus(gadem.sel)

bestmotif <- which(nOccurrences(gadem.sel)==max(nOccurrences(gadem.sel)))

consensus(gadem.seö)[bestmotif]


motifs.sel <- getPWM(gadem.sel)

seqLogo(motifs.sel[[bestmotif]])
seqLogo(reverseComplement(motifs.sel[[bestmotif]]))

#http://www.factorbook.org/mediawiki/index.php/PU.1

#package TFBSTools can be used for further using of PWM
