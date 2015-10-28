library(igraph)


g <- graph(c(1,2, 1,3, 1,4, 2,3, 2,4, 3,4, 4,5, 4,6, 6,7 ),directed=F)
V(g)$name<-c("CBX3","CBX5","CBX1","ARL5A","NMT1","KPNA2","ARL4A")
plot(g)


net <- read.csv(file ="Arabidopsis_Mendosa.sif",sep="\t",header=F)
gd<-graph.edgelist(as.matrix(net[, -2]), directed = T)

E(g)
E(gd)
list.vertex.attributes(g)
get.vertex.attribute(g,"label")
V(g)$label

E(gd)$color<-rep("green",length(E(gd)$type))
E(gd)$color[E(gd)$type=="R"]<-"red"
plot(gd,edge.label=E(gd)$type)

E(g)$weight<-sample(1:4,9,replace=T)
plot(g, edge.width=E(g)$weight)
plot(g, edge.label=E(g)$weight)

sp<-get.shortest.paths(g,from="ARL4A",to="CBX3")
V(g)$name[sp$vpath[[1]]]

get.shortest.paths(gd,from="LUG",to="AP3")
get.shortest.paths(gd,from="AP3",to="LUG")


is.connected(g)
is.connected(graph.empty(10))
is.connected(gd,mode="weak")
is.connected(gd,mode="strong")

par(mfrow=c(2,2))
plot(erdos.renyi.game(200,200,"gnm",directed=F),vertex.label=NA,vertex.size=2,main="ER gnm")
plot(erdos.renyi.game(200,0.01,"gnp",directed=F),vertex.label=NA,vertex.size=2,main="ER gnp")
plot(watts.strogatz.game(1, 200, 1, 0.5),vertex.label=NA,vertex.size=2,main="WS")
plot(barabasi.game(200,directed=F),vertex.label=NA,vertex.size=2,main="BA")
par(mfrow=c(1,1))


############################################################

# read.sif by Leo Lahti <leo.lahti@iki.fi> from netresponse package

read.sif <- function (sif.file, format = "graphNEL", directed = FALSE, header = TRUE, sep = "\t", ...) 
{

    net <- read.csv(file = sif.file, sep = sep, colClasses = "character", header = header, ...)
        
    # Assume form: node1 linktype node2 side.info..
    if ( ncol(net) > 2 ) { 

      # remove NA nodes 
      nas <- apply(net, 1, function (x) {any(is.na(x[c(1,3)]))})
      if (any(nas)) {
        net <- net[!nas, ]
        warning("NAs removed from network node list, ", sum(nas), " edges removed.")
      }
      
      net <- graph.edgelist(as.matrix(net[, -2]), directed = directed)

    } else if ( ncol(net) == 2 ) { # assume form: node1 node2

      # remove NA nodes 
      nas <- apply(net, 1, function (x) {any(is.na(x))})
      if (any(nas)) {
        net <- net[!nas, ]
        warning("NAs removed from network node list, ", sum(nas), " edges removed.")
      }
      
      net <- graph.edgelist(cbind(net[,1],net[,2]), directed = directed)
    }

    if (format == "graphNEL") { net <- igraph.to.graphNEL(net) }
    # if (format == "igraph") { net <- igraph.from.graphNEL(igraph.to.graphNEL(net)) }

    net
}
