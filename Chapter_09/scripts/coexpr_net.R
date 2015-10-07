library(qpgraph)

furin.data<-read.csv("furin_significant_genes.csv",row.names=1)

#data(EcoliOxygen)

pcc <- qpPCC(furin.data)
pcc$Rsign<-pcc$R
pcc$Rsign[pcc$P>0.05]<-NA



nrr.q1 <- qpNrr(as.matrix(furin.data),q=1)
nrr.q3 <- qpNrr(as.matrix(furin.data),q=3)
nrr.q5 <- qpNrr(as.matrix(furin.data),q=5)
nrr.q7 <- qpNrr(as.matrix(furin.data),q=7)

avgnrr <- qpAvgNrr(furin.data)

#density cutoff

par(mfrow=c(3,2))
plot(density(as.numeric(pcc$R),na.rm=T),main="Pearson correlation\nAll")
plot(density(as.numeric(pcc$Rsign),na.rm=T),main="Pearson correlation\nSignificants")
plot(density(as.numeric(nrr.q1),na.rm=T),main="Non-rejection rate\nq=1")
plot(density(as.numeric(nrr.q3),na.rm=T),main="Non-rejection rate\nq=3")
plot(density(as.numeric(nrr.q5),na.rm=T),main="Non-rejection rate\nq=5")
plot(density(as.numeric(nrr.q7),na.rm=T),main="Non-rejection rate\nq=7")
par(mfrow=c(1,1))

par(mfrow=c(2,2))
plot(density(as.numeric(abs(pcc$Rsign)),na.rm=T),main="Pearson correlation\nSignificants")
plot(density(as.numeric(nrr.q1),na.rm=T),main="Non-rejection rate\nq=1",xlim=c(1,0))
plot(density(as.numeric(nrr.q3),na.rm=T),main="Non-rejection rate\nq=3",xlim=c(1,0))
plot(density(as.numeric(nrr.q5),na.rm=T),main="Non-rejection rate\nq=5",xlim=c(1,0))
par(mfrow=c(1,1))



par(mfrow=c(2,2))
qpGraphDensity(nrr.q1, title="q=1", breaks=10)
qpGraphDensity(nrr.q3, title="q=3", breaks=10)
qpGraphDensity(nrr.q5, title="q=5", breaks=10)
qpGraphDensity(nrr.q7, title="q=7", breaks=10)
par(mfrow=c(1,1))


library(graph)

# This is deprecated!
#g <- qpGraph(nrr.q3, threshold=0.1, return.type="graphNEL")
# Workaround based on converting the nrr.q3 object to a matrix, filter
# for the threshold, and convert it to an edgelist (list of edges) and
# coerce a graph object as edgelist

gm <- as.matrix(nrr.q3)
thres <- 0.1
my.nodes <- row.names(gm)
edL <- vector("list", length=length(my.nodes))
names(edL) <- my.nodes
for(i in 1:length(my.nodes)){
  edL[[i]] <- list(edges=names(which(gm[i,]<thres)), weights=gm[i,which(gm[i,]<thres)])
}
g <- graphNEL(nodes=my.nodes, edgeL=edL)

# end of workaround

table(sapply(connComp(g),length))
qpPlotNetwork(g,minimumSizeConnComp=3)

# pcc100 <- qpGraph(as(1-abs(pcc$Rsign),"dspMatrix"),
#                   topPairs=100,
#                   return.type="graphNEL")

gmpc100 <- as.matrix(abs(pcc$Rsign))
gmpc100[gmpc100==1] <- NA
topnum <- 200

pc100vals <- as.vector(gmpc100)
pc100vals <- sort(pc100vals[!is.na(pc100vals)],decreasing = T)

thres <- pc100vals[topnum]

my.nodes.pc100 <- row.names(gmpc100)
edLpc100 <- vector("list", length=length(my.nodes.pc100))
names(edLpc100) <- my.nodes.pc100
for(i in 1:length(my.nodes.pc100)){
  edLpc100[[i]] <- list(edges=names(which(gmpc100[i,]>=thres)), 
                        weights=gmpc100[i,which(gmpc100[i,]>=thres)])
}

pcc100 <- graphNEL(nodes=my.nodes.pc100, edgeL=edLpc100)


table(sapply(connComp(pcc100),length))
qpPlotNetwork(pcc100,minimumSizeConnComp=4)

# Deprecated code:
# qpg100 <- qpGraph(avgnrr,
#                   threshold=NULL,
#                   topPairs=100,
#                   return.type="graphNEL")
gmpg100 <- as.matrix(nrr.q3)
gmpg100[gmpg100==1] <- NA
topnum <- 200

pg100vals <- as.vector(gmpg100)
pg100vals <- sort(pg100vals[!is.na(pg100vals)])

thres <- pg100vals[topnum]

my.nodes.pg100 <- row.names(gmpg100)
edLpg100 <- vector("list", length=length(my.nodes.pg100))
names(edLpg100) <- my.nodes.pg100
for(i in 1:length(my.nodes.pg100)){
  edLpg100[[i]] <- list(edges=names(which(gmpg100[i,]<=thres)), 
                        weights=gmpg100[i,which(gmpg100[i,]<=thres)])
}

qpg100 <- graphNEL(nodes=my.nodes.pg100, edgeL=edLpg100)




table(sapply(connComp(qpg100),length))
qpPlotNetwork(qpg100,minimumSizeConnComp=4)

par(mfrow=c(2,1))
qpPlotNetwork(qpg100,minimumSizeConnComp=6)
qpPlotNetwork(pcc100,minimumSizeConnComp=4)
par(mfrow=c(1,1))


par(mfrow=c(2,1))
qpPlotNetwork(qpg100,vertexSubset="Lhx2", boundary=TRUE)
qpPlotNetwork(pcc100,vertexSubset="Lhx2", boundary=TRUE)
par(mfrow=c(1,1))

