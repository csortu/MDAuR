library(seqinr)
ttn<-read.fasta(file="human_ttn_genomic_dna.fasta")
seq<-ttn[[1]]

table(seq) # We get only nucleotide counts
table(seq)/length(seq) # We get nucleotide frequencies

seq.rho<-rho(seq) # the rho statistics
seq.zsq.base<-zscore(seq, model = "base", exact = F) # using the 'base' model
seq.zsq.codon<-zscore(seq, model = "codon", exact = F) # using the 'codon' model
par(mfrow=c(2,2))
plot(seq.rho,ylim=c(min(seq.rho),max(seq.rho)), lwd = 10, col=col, main="Rho")
plot(seq.rho-1,ylim=c(min(seq.rho)-1,max(seq.rho)-1), lwd = 10, col=col, main="Rho – 1")
plot(seq.zsq.codon,ylim=c(min(seq.zsq.codon),max(seq.zsq.codon)), lwd = 10, col=col, main="Z-score\nCodon model")
plot(seq.zsq.base,ylim=c(min(seq.zsq.base),max(seq.zsq.base)), lwd = 10, col=col, main="Z-score\nBase model")


par(mfrow=c(2,2))

for (win.size in c(1000,5000,10000,50000)){

  gc.perc<-vector()
  k<-1

  for (i in seq(from=1,to=length(seq)-win.size,by=win.size/10)){
    j<-i+win.size-1
    gc.perc[k]<-GC(seq[i:j])
    k<-k+1
  }

  plot(gc.perc*100,type="l", main=win.size,ylab="GC%",xlab="")
}

