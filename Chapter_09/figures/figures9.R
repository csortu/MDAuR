## figure 9.1

library(qpgraph)
furin.data<-read.csv("furin_significant_genes.csv",row.names=1)
pcc <- qpPCC(furin.data)
pcc$Rsign<-pcc$R
pcc$Rsign[pcc$P>0.05]<-NA


nrr.q1 <- qpNrr(as.matrix(furin.data),q=1)
nrr.q3 <- qpNrr(as.matrix(furin.data),q=3)
nrr.q5 <- qpNrr(as.matrix(furin.data),q=5)
nrr.q7 <- qpNrr(as.matrix(furin.data),q=7)

tiff("figure9_1.tif",width=5,height=5,units="in",res=1200)
par(mfrow=c(3,2),mai=c(0.42,0.42,0.42,0.12)) 
plot(density(as.numeric(pcc$R),na.rm=T),
     main="Pearson correlation\nAll",xlab="",ylab="")
mtext("A", side = 3, line = 0, adj = 0, cex = 1)
plot(density(as.numeric(pcc$Rsign),na.rm=T,xlab=""),
     main="Pearson correlation\nSignificants",xlab="",ylab="")
mtext("B", side = 3, line = 0, adj = 0, cex = 1)
plot(density(as.numeric(nrr.q1),na.rm=T),
     main="Non-rejection rate\nq=1",xlab="",ylab="")
mtext("C", side = 3, line = 0, adj = 0, cex = 1)
plot(density(as.numeric(nrr.q3),na.rm=T),
     main="Non-rejection rate\nq=3",xlab="",ylab="")
mtext("D", side = 3, line = 0, adj = 0, cex = 1)
plot(density(as.numeric(nrr.q5),na.rm=T),
     main="Non-rejection rate\nq=5",xlab="",ylab="")
mtext("E", side = 3, line = 0, adj = 0, cex = 1)
plot(density(as.numeric(nrr.q7),na.rm=T),
     main="Non-rejection rate\nq=7",xlab="",ylab="")
mtext("F", side = 3, line = 0, adj = 0, cex = 1)
dev.off()

## figure 9.2

tiff("figure9_2.tif",width=5,height=5,units="in",res=1200)
par(mfrow=c(2,2),mai=c(0.52,0.42,0.52,0.12))
plot(density(as.numeric(abs(pcc$Rsign)),na.rm=T),
     main="Pearson correlation\nSignificants",ylim=c(0,6))
plot(density(as.numeric(nrr.q1),na.rm=T),
     main="Non-rejection rate\nq=1",xlim=c(1,0),ylim=c(0,6))
plot(density(as.numeric(nrr.q3),na.rm=T),
     main="Non-rejection rate\nq=3",xlim=c(1,0),ylim=c(0,6))
plot(density(as.numeric(nrr.q5),na.rm=T),
     main="Non-rejection rate\nq=5",xlim=c(1,0),ylim=c(0,6))
dev.off()

## figure 9.3

tiff("figure9_3.tif",width=5,height=5,units="in",res=1200)
par(mfrow=c(2,2))
qpGraphDensity(nrr.q1, title="q=1", breaks=10)
qpGraphDensity(nrr.q3, title="q=3", breaks=10)
qpGraphDensity(nrr.q5, title="q=5", breaks=10)
qpGraphDensity(nrr.q7, title="q=7", breaks=10)
dev.off()

## figure 9.4
gm <- as.matrix(nrr.q3)
thres <- 0.1
my.nodes <- row.names(gm)
edL <- vector("list", length=length(my.nodes))
names(edL) <- my.nodes
for(i in 1:length(my.nodes)){
  edL[[i]] <- list(edges=names(which(gm[i,]<thres)), weights=gm[i,which(gm[i,]<thres)])
}
library(graph)
g <- graphNEL(nodes=my.nodes, edgeL=edL)





pcc100 <- qpAnyGraph(abs(pcc$Rsign), 
                     threshold=NULL, 
                     topPairs=100, 
                     decreasing=TRUE, 
                     return.type="graphNEL")
qpg100 <- qpGraph(nrr.q3, 
                  topPairs=100, 
                  return.type="graphNEL")


