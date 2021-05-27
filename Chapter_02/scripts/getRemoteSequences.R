library(seqinr)

choosebank()
choosebank(infobank=T)
choosebank("swissprot")
dbg<-get.db.growth()
str(dbg)
head(dbg)
plot(dbg$date,dbg$Nucleotides)

hoxhuman<-query("hoxhuman", 
                "sp=Homo sapiens AND k=homeobox",
                verbose = TRUE)

hoxhuman$nelem # $call $name actual records: $req
names <- getName(hoxhuman$req[1:5])
plot(density(getLength(hoxhuman$req)))
order(getLength(hoxhuman$req)) # human hox proteins from shortest to longest
annot <- getAnnot(hoxhuman$req[2])
seq<-getSequence(hoxhuman$req[2])
table(seq)

# Fetch sequences from NCBI

library(reutils)

accs <- c("NP_002135",
          "NP_001074229",
          "NP_001075041.1",
          "XP_548172.3",
          "NP_001192883.1",
          "NP_032292.3",
          "XP_220896.4")

fastafile <- 'new_sequences.fst'

mseq <- efetch(accs, db = "protein",
               rettype = 'fasta',
               retmode = 'text',
               outfile = fastafile)

my.prots <- read.fasta(fastafile,seqtype = 'AA')
my.annot <- unlist(getAnnot(my.prots))
my.orgs <- gsub("\\]","",gsub("^.*\\[","",unlist(getAnnot(my.prots))))

AAstat(unlist(getSequence(my.prots[1])))

pdf("class5_aastats.pdf")
par(mfrow=c(3,3))
for(i in 1:length(my.prots)){
  AAstat(unlist(getSequence(my.prots[1])))
  title(my.orgs[i])
}
dev.off()
