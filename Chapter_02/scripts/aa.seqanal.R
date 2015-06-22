library(seqinr)
rhod<-read.fasta("human_rhodopsin_protein.fasta",seqtype = "AA")
rhodseq<-rhod[[1]]
AAstat(rhodseq)

hpat<-read.table("aa_hydropathy.txt")
hpatseq<-hpat[rhodseq[1:length(rhodseq)],2]

par(mfrow=c(2,2))

for (win.size in c(2,4,8,16)){

  hpat.m<-vector()
  k<-1

  for (i in seq(from=1,to=length(hpatseq)-win.size)){
    j<-i+win.size-1
    hpat.m[k]<-mean(hpatseq[i:j])
    k<-k+1
  }

  plot(hpat.m,type="l", main=win.size,ylab="Hydropathy",xlab="Sequence")
}

