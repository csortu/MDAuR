library(igraph)
gw<-read.graph("PredictedInteractions1000000.ncol",format="ncol")
V(gw)$cluster<-clusters(gw)$membership
g.human<-simplify(delete.vertices(gw,V(gw)$cluster!=2))
#g.yeast<-read.graph("yeast.graphml",format="graphml")

plot(g.human)
plot(g.yeast)

plot(g.human,vertex.label=NA, vertex.size=sample(1:5, vcount(g.human), replace=T))

list.edge.attributes(g.human)
E(g.human)$lty<-1
E(g.human)$lty[log(E(g.human)$weight)<mean(log(E(g.human)$weight))]<-3
plot(g.human,vertex.label=NA,vertex.size=5)

par(mfrow=c(2,2))
plot(g.human,vertex.label=NA,vertex.size=5,layout=layout.random)
plot(g.human,vertex.label=NA,vertex.size=5,layout=layout.circle)
plot(g.human,vertex.label=NA,vertex.size=5,layout=layout.kamada.kawai)
plot(g.human,vertex.label=NA,vertex.size=5,layout=layout.reingold.tilford)
par(mfrow=c(1,1))

g.hu.coords<-layout.kamada.kawai(g.human)
plot(g.human,vertex.label=NA,vertex.size=5,layout=g.hu.coords)

V(g.human)$x<-g.hu.coords[,1]
V(g.human)$y<-g.hu.coords[,2]
par(mfrow=c(2,2))
plot(g.human,vertex.label=NA, vertex.size=sample(1:5, vcount(g.human), replace=T))
plot(g.human,vertex.label=NA, vertex.color=sample(c("red","green","blue"), vcount(g.human), replace=T),edge.lty=1)
plot(g.human,vertex.label=NA, edge.color=sample(c("red","green","blue"), ecount(g.human), replace=T),edge.lty=1,edge.width=3,vertex.size=4)
plot(g.human,vertex.label=NA, vertex.size=sample(1:5, vcount(g.human), replace=T),edge.label=E(g.human)$weight,edge.lty=1)
par(mfrow=c(1,1))

tkplot(g.human,vertex.label=NA, edge.color=sample(c("red","green","blue"), ecount(g.human), replace=T),edge.lty=1,edge.width=3,vertex.size=4)

plotindex<-tkplot(g.human,vertex.label=NA, edge.color=sample(c("red", "green", "blue"), ecount(g.human),  replace=T), edge.lty=1,edge.width=3,vertex.size=4)
tkplot.getcoords(plotindex)

library(igraph)
library(RCytoscape)

net <- read.csv(file ="Arabidopsis_Mendosa.sif",sep="\t",header=F)
g.arabid<-graph.edgelist(as.matrix(net[, -2]), directed = T)
E(g.arabid)$type<-as.character(net[,2])

g<-igraph.to.graphNEL(g.arabid)
g<-initEdgeAttribute (g,'type','char',default.value='A')
g<-initEdgeAttribute (g,'weight','numeric',default.value=1)
g<-initNodeAttribute (g,'label','char',default.value='Gene')

cw <- new.CytoscapeWindow ('Arabidopsis GRN', graph=g)
displayGraph (cw)
layoutNetwork (cw, 'jgraph-spring')

