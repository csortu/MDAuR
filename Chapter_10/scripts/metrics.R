library(igraph)
net <- read.csv(file ="Arabidopsis_Mendosa.sif",sep="\t",header=F)
g.arabid<-graph.edgelist(as.matrix(net[, -2]), directed = T)
E(g.arabid)$type<-as.character(net[,2])

gw<-read.graph("PredictedInteractions1000000.ncol",format="ncol")
V(gw)$cluster<-clusters(gw)$membership
g.human<-delete.vertices(gw,V(gw)$cluster!=2)
g.human<-simplify(g.human)

g.yeast<-read.graph("yeast.graphml",format="graphml")

gr.erdos<-erdos.renyi.game(200,200,"gnm",directed=F)
V(gr.erdos)$cluster<-clusters(gr.erdos)$membership
gr.erdos<-delete.vertices(gr.erdos,V(gr.erdos)$cluster!=1)

gr.watts<-watts.strogatz.game(1, 200, 1, 0.5)
V(gr.watts)$cluster<-clusters(gr.watts)$membership
gr.watts<-delete.vertices(gr.watts,V(gr.watts)$cluster!=1)

gr.barabasi<-barabasi.game(200,directed=F)

par(mfrow=c(2,3))
plot(g.arabid, main="Arabidopsis GNR")
plot(g.human,main="Human PPI",vertex.size=2)
plot(g.yeast,main="Yeast PPI",vertex.size=2)
plot(gr.erdos,main="ER random graph",vertex.size=2, vertex.label=NA)
plot(gr.watts,main="WS random graph",vertex.size=2, vertex.label=NA)
plot(gr.barabasi,main="BA random graph",vertex.size=2, vertex.label=NA)
par(mfrow=c(1,1))


shortest.paths(g.arabid,"LUG","AP3")
get.shortest.paths(g.arabid,"LUG","AP3")

shortest.paths(g.human,"ZAP70","JUN")
get.shortest.paths(g.human,"ZAP70","JUN")

shortest.paths(g.human,"ZAP70","JUN", weights=NA)
get.shortest.paths(g.human,"ZAP70","JUN", weights=NA)

shortest.paths(g.arabid,"LUG")
shortest.paths(g.arabid)

par(mfrow=c(2,3))
plot(density(shortest.paths(g.arabid)), main="Arabidopsis GNR")
plot(density(shortest.paths(g.human,weights=NA)), main="Human PPI")
plot(density(shortest.paths(g.yeast)), main="Yeast PPI")
plot(density(shortest.paths(gr.erdos)), main="ER random graph")
plot(density(shortest.paths(gr.watts)), main="WS random graph")
plot(density(shortest.paths(gr.barabasi)), main="BA random graph")
par(mfrow=c(1,1))

average.path.length(g.human)
average.path.length(gr.erdos)
average.path.length(gr.watts)
average.path.length(gr.barabasi)

global.efficiency<- function(graph){

  if (!is.igraph(graph)){stop("The parameter is not an igraph object!")}

  N<-vcount(graph)
  dij<-shortest.paths(graph)
  dij[dij==0]<-Inf

  out<-(1/(N*(N-1))) * sum(1/dij)
  return(out)
}

barplot(c(global.efficiency(g.arabid),global.efficiency(g.human),global.efficiency(g.yeast),global.efficiency(gr.erdos),global.efficiency(gr.watts),global.efficiency(gr.barabasi)),names.arg=c("Arabidopsis","Human","Yeast","ER","WS","BA"),ylab="Global efficiency")

degree(g.arabid,"LFY")
degree(g.arabid,"LFY",mode="in")
degree(g.arabid,"LFY",mode="out")

degree(g.arabid)
degree(g.arabid,mode="in")
degree(g.arabid,mode="out")

degree.distribution(g.arabid)

par(mfrow=c(2,3))
plot(degree.distribution(g.arabid), type="l", ylab="", xlab="Degree", main="Arabidopsis")
plot(degree.distribution(g.human), type="l", ylab="", xlab="Degree", main="Human")
plot(degree.distribution(g.yeast), type="l", ylab="", xlab="Degree", main="Yeast")
plot(degree.distribution(gr.erdos), type="l", ylab="", xlab="Degree", main="ER")
plot(degree.distribution(gr.watts), type="l", ylab="", xlab="Degree", main="WS")
plot(degree.distribution(gr.barabasi), type="l", ylab="", xlab="Degree", main="BA")
par(mfrow=c(1,1))


vulnerability<-function(graph){
  if (!is.igraph(graph)){stop("The parameter is not an igraph object!")}
  E<-global.efficiency(graph)
  vuln<-vector()
  for (node in V(graph)){
    gr.del<-delete.vertices(graph,node)
    Ei<-global.efficiency(gr.del)
    Vi<-(E-Ei)/E
    vuln<-c(vuln,Vi)
  }
  return(vuln)
}

V(g.arabid)$vuln<-vulnerability(g.arabid)
max(V(g.arabid)$vuln)

g.arabid$vulnarability<-max(V(g.arabid)$vuln)

library(gplots)
V(g.arabid)$color<-rev(redgreen(100))[cut(V(g.arabid)$vuln, 100)]
plot(g.arabid)

V(g.arabid)$vuln<-vulnerability(g.arabid)
V(g.arabid)$color<-rev(redgreen(100))[cut(V(g.arabid)$vuln, 100)]

V(g.human)$vuln<-vulnerability(g.human)
V(g.human)$color<-rev(redgreen(100))[cut(V(g.human)$vuln, 100)]

V(gr.erdos)$vuln<-vulnerability(gr.erdos)
V(gr.erdos)$color<-rev(redgreen(100))[cut(V(gr.erdos)$vuln, 100)]

V(gr.barabasi)$vuln<-vulnerability(gr.barabasi)
V(gr.barabasi)$color<-rev(redgreen(100))[cut(V(gr.barabasi)$vuln, 100)]


par(mfrow=c(2,2))
plot(g.arabid, main="Arabidopsis GNR")
plot(g.human,main="Human PPI",vertex.size=5, vertex.label=NA)
plot(gr.erdos,main="ER random graph",vertex.size=5, vertex.label=NA)
plot(gr.barabasi,main="BA random graph",vertex.size=5, vertex.label=NA)
par(mfrow=c(1,1))

transitivity(g.human,type="global")
V(g.human)$cc<-transitivity(g.human,type="local")
V(g.human)$color<-rev(redgreen(100))[cut(V(g.human)$cc, 100)]
plot(g.human,main="Human PPI",vertex.size=5, vertex.label=NA)

fc.b<-fastgreedy.community(gr.barabasi)
plot(fc.b,gr.barabasi,vertex.label=NA, vertex.size=5)

fc.h<-fastgreedy.community(g.human)
plot(fc.h,g.human,vertex.label=NA, vertex.size=5)

wc.h<-walktrap.community(g.human)
plot(wc.h,g.human,vertex.label=NA, vertex.size=5)

sc.h<-spinglass.community(g.human)
plot(sc.h,g.human,vertex.label=NA, vertex.size=5)

V(g.human)[V(g.human)$sc==V(g.human)$sc[V(g.human)$name=="JUN"]]

