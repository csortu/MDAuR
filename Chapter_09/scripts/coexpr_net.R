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


g <- qpGraph(nrr.q1, threshold=0.001, return.type="graphNEL")
library(graph)
table(sapply(connComp(g),length))
qpPlotNetwork(g,minimumSizeConnComp=3)

pcc100 <- qpAnyGraph(abs(pcc$Rsign), threshold=NULL,topPairs=100, decreasing=TRUE, return.type="graphNEL")
table(sapply(connComp(pcc100),length))
qpPlotNetwork(pcc100,minimumSizeConnComp=4)

qpg100 <- qpGraph(avgnrr, threshold=NULL, topPairs=100, return.type="graphNEL")
table(sapply(connComp(qpg100),length))
qpPlotNetwork(qpg100,minimumSizeConnComp=6)

par(mfrow=c(2,1))
qpPlotNetwork(qpg100,minimumSizeConnComp=6)
qpPlotNetwork(pcc100,minimumSizeConnComp=4)
par(mfrow=c(1,1))


par(mfrow=c(2,1))
qpPlotNetwork(g,vertexSubset="Lhx2", boundary=TRUE)
qpPlotNetwork(qpg100,vertexSubset="Lhx2", boundary=TRUE)
par(mfrow=c(1,1))

