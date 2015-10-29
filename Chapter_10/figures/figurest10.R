library(igraph)

## figure 10.1

g0 <- graph(c(1,2,1,3,1,4,2,3),directed=T)
V(g0)$color <- 'white'
V(g0)$label.color <- 'black'

g.layout <- matrix(c(0.40293486,-0.03631870,
                    -0.14276948,0.08694801,
                    0.60070463,-0.73359732,
                    0.08169628,0.40864797),
                   ncol=2,byrow=T)

tiff("figure10_1.tif",width=5,height=5,units="in",res=1200)
plot(g0,layout=g.layout,
     vertex.size=30,
     vertex.label=c("node","","",""),
     vertex.label.cex=1.6,
     vertex.shape=c("none","circle","circle","circle"),
     edge.color="black",
     edge.label=c("15","","","4"),
     edge.label.color="black",edge.width=c(7,2,2,4),
     edge.label.cex=2,
     edge.label.x=c(-0.2,NA,NA,0.15),
     edge.label.y=c(0.45,NA,NA,-0.15),
     edge.arrow.mode=c(2,2,0,0),
     edge.arrow.width=2)
dev.off()

## figure 10.1

g <- graph(c(1,2, 1,3, 1,4, 2,3, 2,4, 3,4, 4,5, 4,6, 6,7 ), directed=F)

g.layout <- matrix(c(-0.02674107,2,
                     0.10762234,0.89263332,
                     0.56276658,1.99632320,
                     3.02243711,0.63438828,
                     1.39691511,1.37579132,
                     0.41504407,0,
                     0.30007442,0.21528329
                     ),nrow=7,ncol=2,byrow = F)

V(g)$color <- 'white'
V(g)$label.color <- 'black'

tiff("figure10_2.tif",width=5,height=5,units="in",res=1200)
par(mai=c(0,0,0,0))
V(g)$name<-c("CBX3","CBX5","CBX1","ARL5A","NMT1","KPNA2","ARL4A")
plot(g,layout=g.layout,vertex.shape='none')
dev.off()

## figure 10.3

net <- read.csv(file ="Arabidopsis_Mendosa.sif",sep="\t",header=F)
gd<-graph.edgelist(as.matrix(net[, -2]), directed = T)
V(gd)$color <- 'white'
V(gd)$label.color <- 'black'

gd.layout <- layout.kamada.kawai(gd)


tiff("figure10_3.tif",width=5,height=5,units="in",res=1200)
par(mai=c(0,0,0,0))
plot(gd,layout=gd.layout,
     vertex.shape='none',
     edge.arrow.width=0.5,
     edge.color="gray")
dev.off()

## figure 10.4

E(gd)$type<-as.character(net[,2])
E(gd)$color<-rep("gray",length(E(gd)$type))
E(gd)$color[E(gd)$type=="R"]<-"black"

tiff("figure10_4.tif",width=5,height=5,units="in",res=1200)
par(mai=c(0,0,0,0))
plot(gd,edge.label=E(gd)$type,
     layout=gd.layout,
     vertex.shape='none',
     edge.arrow.width=0.5,
     edge.label.color="black")
dev.off()

## figure 10.5

tiff("figure10_5.tif",width=5,height=5,units="in",res=1200)
par(mfrow=c(2,2),mai=c(0,0,0.4,0))
plot(erdos.renyi.game(200,200,"gnm",directed=F),
     vertex.label=NA,vertex.size=2,
     main="Erdős-Rényi (gnm)")
plot(erdos.renyi.game(200,0.01,"gnp",directed=F), 
     vertex.label=NA, vertex.size=2, 
     main="Erdős-Rényi (gnp)")
plot(watts.strogatz.game(1, 200, 1, 0.5), 
     vertex.label=NA, vertex.size=2, 
     main="Watts-Strogatz")
plot(barabasi.game(200,directed=F), 
     vertex.label=NA, vertex.size=2, 
     main="Barabási-Albert")
dev.off()

## figure 10.6

gw<-read.graph("PredictedInteractions1000000.ncol",format="ncol")

tiff("figure10_6.tif",width=5,height=5,units="in",res=1200)
par(mfrow=c(1,2),mai=c(0,0,0.4,0))
plot(gw,     
     vertex.size=3,
     vertex.color="white",
     vertex.label=NA,
     edge.width=1.5,
     edge.color="black")
mtext("A", side = 3, line = -2, adj = 0.01, cex = 1.3)
V(gw)$cluster<-clusters(gw)$membership
gws<-delete.vertices(gw,V(gw)$cluster!=2)
plot(gws,
     vertex.size=3,
     vertex.color="white",
     vertex.label=NA,
     edge.width=1.5,
     edge.color="black")
mtext("B", side = 3, line = -2, adj = 0.01, cex = 1.3)
dev.off()

## figure 10.7

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


tiff("figure10_7.tif",width=5,height=5,units="in",res=1200)
par(mfrow=c(2,3),mai=c(0,0,0.4,0))
plot(g.arabid, 
     main="Arabidopsis GRN",
     vertex.label.color="black",
     vertex.shape='none',
     edge.arrow.width=0.4,
     edge.arrow.size=0.5,
     edge.color="gray")
plot(g.human,main="Human PPI",
     vertex.label.color="black",
     vertex.label.cex=0.5,
     vertex.shape="none",
     edge.width=1.5,
     edge.color="black")
plot(g.yeast,main="Yeast PPI",
     vertex.size=2,
     vertex.color="black",
     vertex.label=NA)
plot(gr.erdos,main="Erdos random graph", 
     vertex.size=2, 
     vertex.color="black",
     vertex.label=NA)
plot(gr.watts,main="WS random graph", 
     vertex.size=2, 
     vertex.color="black",
     vertex.label=NA)
plot(gr.barabasi,main="BA random graph", 
     vertex.size=2, 
     vertex.color="black",
     vertex.label=NA)
dev.off()

## figure 10.8

global.efficiency<- function(graph){
  if (!is.igraph(graph)){stop("The parameter is not an igraph object!")}
  N<-vcount(graph)
  dij<-shortest.paths(graph)
  dij[dij==0]<-Inf
  out<-(1/(N*(N-1)))*sum(1/dij)
  return(out)
}

geff.arabid <- global.efficiency(g.arabid)
geff.human <- global.efficiency(g.human)
geff.yeast <- global.efficiency(g.yeast)
geff.erdos <- global.efficiency(gr.erdos)
geff.watts <- global.efficiency(gr.watts)
geff.barabasi <- global.efficiency(gr.barabasi)

tiff("figure10_8.tif",width=5,height=5,units="in",res=1200)
par(mfrow=c(1,1),mai=c(0.5,0.8,0.4,0))
barplot(c(geff.arabid, 
          geff.human, 
          geff.yeast, 
          geff.erdos, 
          geff.watts, 
          geff.barabasi), 
        names.arg=c("Arab.","Human","Yeast","ER","WS","BA"), 
        ylab="Global efficiency")
dev.off()

## figure 10.9

tiff("figure10_9.tif",width=5,height=5,units="in",res=1200)
par(mfrow=c(2,3),mai=c(0.3,0.25,0.4,0.05))
plot(degree.distribution(g.arabid), 
     type="l", ylab="", main="Arabidopsis")
plot(degree.distribution(g.human), 
     type="l", ylab="", main="Human")
plot(degree.distribution(g.yeast), 
     type="l", ylab="", main="Yeast")
plot(degree.distribution(gr.erdos), 
     type="l", ylab="", main="ER")
plot(degree.distribution(gr.watts), 
     type="l", ylab="", main="WS")
plot(degree.distribution(gr.barabasi), 
     type="l", ylab="", main="BA")
dev.off()

## figure 10.10

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
V(g.arabid)$color <- rev(gray(1:100/100))[cut(V(g.arabid)$vuln, 100)]
V(g.arabid)$label.color <- gray(1:100/100)[cut(V(g.arabid)$vuln, 100)]

tiff("figure10_10.tif",width=5,height=5,units="in",res=1200)
par(mai=c(0,0,0,0))
plot(g.arabid,
     vertex.shape='rectangle',
     vertex.size=30,vertex.size2=15,
     edge.arrow.width=0.5,
     edge.color="black")
dev.off()

## figure 10.11

gw<-read.graph("PredictedInteractions1000000.ncol",format="ncol")
V(gw)$cluster<-clusters(gw)$membership
g.human<-delete.vertices(gw,V(gw)$cluster!=2)
g.human<-simplify(g.human)

V(g.human)$cc<-transitivity(g.human,type="local")
#V(g.human)$color<-rev(redgreen(100))[cut(V(g.human)$cc, 100)]
V(g.human)$color <- "white"
V(g.human)[which(V(g.human)$cc>0.6)]$color <- "black"

# tiff("figure10_11.tif",width=5,height=5,units="in",res=1200)
# par(mai=c(0,0,0,0))
# plot(g.human,vertex.size=5, vertex.label=NA)
# dev.off()


g.human$layout <- layout.fruchterman.reingold(g.human)
V(g.human)$color <- "white"

fc.b<-fastgreedy.community(g.human)
wc.h<-walktrap.community(g.human)
sc.h<-spinglass.community(g.human)

tiff("figure10_11_color.tif",width=5,height=5,units="in",res=1200)
par(mfrow=c(2,2),mai=c(0,0,0,0))
plot(g.human,vertex.label=NA, vertex.size=5)
mtext("A", side = 3, line = -2, adj = 0.01, cex = 1.3)
plot(fc.b,g.human,vertex.label=NA, vertex.size=5)
mtext("B", side = 3, line = -2, adj = 0.01, cex = 1.3)
plot(wc.h,g.human,vertex.label=NA, vertex.size=5)
mtext("C", side = 3, line = -2, adj = 0.01, cex = 1.3)
plot(sc.h,g.human,vertex.label=NA, vertex.size=5)
mtext("D", side = 3, line = -2, adj = 0.01, cex = 1.3)
dev.off()

## figure 10.12

gw<-read.graph("PredictedInteractions1000000.ncol",format="ncol")
V(gw)$cluster<-clusters(gw)$membership
g.human<-simplify(delete.vertices(gw,V(gw)$cluster!=2))
V(g.human)$color <- "grey"

g.human$layout <- layout.fruchterman.reingold(g.human)

tiff("figure10_12.tif",width=5,height=5,units="in",res=1200)
par(mfrow=c(1,1),mai=c(0,0,0,0))
plot(g.human,vertex.label=NA, vertex.size=sample(1:5, vcount(g.human), replace=T))
dev.off()

## figure 10.13

E(g.human)$lty<-1
E(g.human)$lty[log(E(g.human)$weight)<mean(log(E(g.human)$weight))]<-3

tiff("figure10_13.tif",width=5,height=5,units="in",res=1200)
par(mfrow=c(1,1),mai=c(0,0,0,0))
plot(g.human,vertex.label=NA,vertex.size=5,edge.color="black",edge.width=3)
dev.off()

## figure 10.14

tiff("figure10_14.tif",width=5,height=5,units="in",res=1200)
par(mfrow=c(2,2),mai=c(0,0,0,0))
plot(g.human,vertex.label=NA,vertex.size=5,layout=layout.random)
mtext("A", side = 3, line = -2, adj = 0.01, cex = 1.3)
plot(g.human,vertex.label=NA,vertex.size=5,layout=layout.circle)
mtext("B", side = 3, line = -2, adj = 0.01, cex = 1.3)
plot(g.human,vertex.label=NA,vertex.size=5,layout=layout.kamada.kawai)
mtext("C", side = 3, line = -2, adj = 0.01, cex = 1.3)
plot(g.human,vertex.label=NA,vertex.size=5,layout=layout.reingold.tilford)
mtext("D", side = 3, line = -2, adj = 0.01, cex = 1.3)
dev.off()
