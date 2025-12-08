library(seqinr)

carps<-read.alignment("vertebrate_CARPs.aln",format="clustal", oldclustal = TRUE)
# carps<-read.alignment("vertebrate_CARPs.aln",format="clustal") # Old seqinr version

carps$nb # The number of sequences in the alignment
carps$nam # The names of the sequences in a vector
carps$seq # A list of sequences in the same order as the names in the carps$nam vector

dali<-dist.alignment(carps)
plot(density((dali)))

carps.res<-kaks(carps)
kaks<-as.matrix(carps.res$ka/carps.res$ks)
density(kaks) # Will throw error

length(kaks)
sum(is.na(kaks))

density(kaks,na.rm=T)
plot(density(kaks,na.rm=T))
