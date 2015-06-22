library(ShortRead)
projectDir <- '/data1/Omixon/target_brca_example'
#setwd(fastqDir)
#fastqFiles <- dir(fastqDir, full=TRUE)
fastqFiles <- list.files(path=projectDir,pattern="*.fastq",full.names=T)

fq <- readFastq(fastqFiles[1])
fq

class(fq) # we will use ShortReadQ class to handle the data


sread(fq) # sequences
#head(sread(fq), 3)
quality(fq)  # quality strings
id(fq) # ids
detail(fq)  # all info at once

width(fq) # length of all reads


# if the file too big use only a sample

sampler <- FastqSampler(fastqFiles[1], 10000) # This will result 10000 random reads

smpl <- yield(sampler)
detail(smpl)
table(as.vector(unlist(sread(smpl))))

# Some basic diagnostics

sread(fq)[1]
as.character(sread(fq)[1])
as.vector(sread(fq)[[1]])
as.vector(sread(fq)[[1]])

table(as.vector(sread(fq)[[1]]))
table(as.vector(quality(fq)[[1]]))

table(as.vector(unlist(sread(fq))))

############################################################
# Quality control

library(ShortRead)
fastqDir <- '/data1/Omixon/target_brca_example'
fastqFiles <- list.files(path=fastqDir,pattern="*.fastq",full.names=T)
fq <- readFastq(fastqFiles[1])


alphabetFrequency(sread(fq), baseOnly=TRUE, collapse=TRUE)
head(alphabetFrequency(sread(fq)))
matplot(t(alphabetByCycle(sread(fq))[c("A","C","G","T"),]),type="l")

colSums(dinucleotideFrequency(sread(fq)))
colSums(trinucleotideFrequency(sread(fq)))

## Get the less and most represented 6-mers:
f6 <- oligonucleotideFrequency(sread(fq), 6)
f6[f6 == min(f6)]
f6[f6 == max(f6)]

# dinucleotide bias analysis

mono <- alphabetFrequency(sread(fq), baseOnly=TRUE, collapse=TRUE)
di <- colSums(dinucleotideFrequency(sread(fq)))

mono.freq <- mono[-5] / sum(mono[-5])
di.exp.freq <- as.vector(outer(mono.freq, mono.freq))
names(di.exp.freq) <- names(di)
di.exp.count <- di.exp.freq * sum(di)
mode(di.exp.count) <- "integer" # round to integer


plot(di.exp.count,di,xlab="Expected count",ylab="Observed count",t="n")
text(di.exp.count,di,labels=names(di))
abline(0,1,lty=2)

# automated quality analysis

# single fastq file
sampler <- FastqSampler(fastqFiles[1], 10000) # This will result 10000 random reads
qa1 <- qa(yield(sampler),basename(fastqFiles[1]))
qa1
qa1[["readCounts"]]
qa1[["baseCalls"]]
head(qa1[["readQualityScore"]])
qa1[["baseQuality"]]
plot(qa1[["baseQuality"]]$count,t="l")

rurl1 <- report(qa1,dest="quality_example1")
browseURL(rurl1)

# all fastq files in directory

qas0 <- Map(function(fl, nm) {
  fq <- FastqSampler(fl)
  qa(yield(fq), nm)
}, fastqFiles, sub(".fastq", "", basename(fastqFiles)))
qas <- do.call(rbind, qas0)
rpt <- report(qas, dest="quality_example")

#browseURL(rpt)

## t.cyc <- Map(function(i,j){
##   print(paste(i,j))
## },rep(1:5,each=2),6:7)



# Quality filtering of fastq files

fq <- readFastq(fastqFiles[1])

# filter for info in the ID line: strand info and fragment ID

fq.plus <- fq[idFilter("Strand \\+")(fq)]
fq3000 <- fq[idFilter("Frag_3[0123456789]{3}")(fq)]
id(fq3000[c(1,1001)])

# filter for unresolved bases: Ns

fq.n <- fq[nFilter()(fq)]
fq
fq.n

# trim low quality tails

#fq.trimmed <- trimTailw(fq.n, 2, ">", 2)
fq.trimmed <- trimTails(fq.n, 5, "A", successive=T)
fq.trimmed
plot(density(width(fq.trimmed)))
 
# filter for length

fq.filtered <- fq.trimmed[width(fq.trimmed) > 80]
plot(density(width(fq.filtered)))

# write out the filtered sequences into a new file

writeFastq(fq.filtered,"brca.example.illumina.filtered.0.1.fastq",compress=F)

# Now do the same for very large files in chunks using a file stream

trim.file <- function(fl, destination=sprintf("%s_filtered.fastq", fl)){

  stream <- open(FastqStreamer(fl))
  on.exit(close(stream))

  repeat {

    fq <- yield(stream)
    if (length(fq) == 0){break}
    fq <- fq[nFilter()(fq)]
    fq <- trimTails(fq, 5, "A", successive=T)
    fq <- fq[width(fq) > 80]
    writeFastq(fq, destination, "a",compress=F )
  }
}

trim.file(fastqFiles[1])


## > myFilterAndTrim <-
## + function(fl, destination=sprintf("%s_subset", fl))
## + {
## + ## open input stream
## + stream <- open(FastqStreamer(fl))
## + on.exit(close(stream))
## +
## + repeat {
## + ## input chunk
## + fq <- yield(stream)
## + if (length(fq) == 0)
## + break
## +An Introduction to ShortRead 5
## + ## trim and filter, e.g., reads cannot contain ✬N✬...
## + fq <- fq[nFilter()(fq)] # see ?srFilter for pre-defined filters
## + ## trim as soon as 2 of 5 nucleotides has quality encoding less
## + ## than "4" (phred score 20)
## + fq <- trimTailw(fq, 2, "4", 2)
## + ## drop reads that are less than 36nt
## + fq <- fq[width(fq) < 36] # should be  >  36 !!!!
## +
## + ## append to destination
## + writeFastq(fq, destination, "a")
## + }
## + }
