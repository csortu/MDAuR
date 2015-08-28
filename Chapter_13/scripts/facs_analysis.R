############################################################
# Simple reading
############################################################
library(flowCore)
fcs<-read.FCS("Migr 1_5h stim IL2.fcs")
fcs
des<-description(fcs)
names(des)
des$FSCversion
des$"$DATE"
des$"$P3N"
des$"$P3V"

summary(fcs)
head(fcs)
exprs(fcs)

library(flowViz)
plot(fcs) # plot all channel pairs
plot(fcs,c("FSC-A","SSC-A")) # Plot a selected pair of channels
plot(fcs,"FSC-A") # Histogram of a single channel

############################################################
# Reading all FCS files
############################################################

library(flowCore)
library(flowViz)

fcspath <- "/home/ortutay/Dropbox/Misc/Zsuzsi/FACS_example/testing"
setwd(fcspath)

dat <- read.flowSet(list.files(fcspath,pattern = ".fcs", full = T))

############################################################
# Setting up phenoData
############################################################

pData(dat)<-data.frame(name=pData(dat)[, 1],infection=NA,stimulation=NA,row.names=row.names(pData(dat)),desc=NA)
varMetadata(phenoData(dat))$labelDescription <- c("File","Infection","Stimulation","Description")
#pData(dat)$infection<-c(rep(NA,4),rep("Cre1",3),rep("Cre2",3),rep("Cre3",3),rep("Migr1",3),rep("Migr2",3),rep("Migr3",3),NA)
pData(dat)$infection<-c(rep(NA,4),rep("Cre",9),rep("Migr",9),NA)
pData(dat)$stimulation<- c(rep(NA,4),rep(c("24h","05h","00h"),6),NA)
pData(dat)$desc<-fsApply(dat,function(f) f@description $`TUBE NAME`)


############################################################
# Compensation matrix
############################################################


namePattern<-"FSC-A|SSC-A|GFP|PE|PerCP-Cy5-5"
spill<-spillover(x=dat[1:4],unstained=4,stain_match="ordered",patt=namePattern,fsc="FSC-A",ssc="SSC-A",method="mean",useNormFilt=T) # We use the mean instead of the default median!

# Check how the calculated spill matrix correlates with the built in spill mtx:

cor.test(as.vector(des$"SPILL"),as.vector(spill))

fcs.compensated<-compensate(fcs,spill)

my.dat<-dat[5:22]

my.dat.compensated<-compensate(my.dat,spill)

xyplot(`SSC-A` ~ `FSC-A` , data = my.dat, subset = infection == "Cre")

############################################################
# Quality control 
############################################################

library(flowQ)

dest <- file.path(tempdir(), "flowQ")

qp.cell <- qaProcess.cellnumber(my.dat, outdir=dest, cFactor=0.75)
qp.margin <- qaProcess.marginevents(my.dat, outdir=dest,pdf=FALSE)
qp.timeline <- qaProcess.timeline(my.dat, channel="PerCP-Cy5-5-A", outdir=dest, cutoff=1)
qp.timeflow <- qaProcess.timeflow(my.dat, outdir=dest, cutoff=2)

url <- writeQAReport(my.dat, processes=list(qp.cell, qp.margin, qp.timeline, qp.timeflow), outdir=dest)

#browseURL(url)


############################################################
# Handle margin events
############################################################

# max values in FSC and SSC channels

max.fsc<-max(exprs(my.dat[[1]])[,"FSC-A"])
max.ssc<-max(exprs(my.dat[[1]])[,"SSC-A"])

# Margin events with max values

margin.cells<-as.logical((exprs(my.dat[[1]])[,"SSC-A"]==max.ssc) + (exprs(my.dat[[1]])[,"FSC-A"]==max.fsc))

# We can plot them

plot(my.dat[[1]],c("FSC-A","SSC-A"),smooth=F)
points(exprs(my.dat[[1]])[margin.cells,c("FSC-A","SSC-A")],col="red")

# We can see what fraction is that in the total count

sum(margin.cells)/nrow(my.dat[[1]])

# remove them from the channel
my.dat[[1]]<-my.dat[[1]][!margin.cells] # negate logical matrix with ! operator

# do for all samples in the set !!Do not run this more than once!!

for (i in 1:18){
  max.fsc<-max(exprs(my.dat[[i]])[,"FSC-A"])
  max.ssc<-max(exprs(my.dat[[i]])[,"SSC-A"])
  margin.cells<-as.logical((exprs(my.dat[[i]])[,"SSC-A"]==max.ssc) + (exprs(my.dat[[i]])[,"FSC-A"]==max.fsc))
  print(pData(my.dat)$name[i])
  print(paste("Margin cell ratio",sum(margin.cells)/nrow(my.dat[[i]])))
  my.dat[[i]]<-(my.dat[[i]])[!margin.cells]
}


############################################################
# Create workflow and add compensation
############################################################

wf <- workFlow(my.dat)
add(wf, compensation(spill, compensationId = "Compensated Samples"))

xyplot(`SSC-A` ~ `FSC-A` , data = wf[["Compensated Samples"]], subset = infection == "Cre")


############################################################
# Transformation, using common Logicle
############################################################

tr<-estimateLogicle(flowFrame(fsApply(my.dat, function(x) exprs(x))),colnames(my.dat)[1:5])
identifier(tr) <- "commonLogicle"

add(wf, tr, parent = "Compensated Samples")

plot(Data(wf[["commonLogicle"]])[[1]])
xyplot(`SSC-A` ~ `FSC-A` , data = wf[["commonLogicle"]], subset = infection == "Cre",smooth=F)

library(flowStats)
densityplot(wf[["base view"]])
X11()
densityplot(wf[["commonLogicle"]])

plot(Data(wf[["base view"]])[[1]],c("FSC-A","SSC-A"))
plot(Data(wf[["commonLogicle"]])[[1]],c("FSC-A","SSC-A"),xlim=c(2.8,4.6),ylim=c(2.8,4.6))


############################################################
# Filtering
############################################################


# Simple rectangle gate

rg <- rectangleGate(filterId="Rectangle","FSC-A"=c(3.6, 4.2), "SSC-A"=c(3.2, 4))
add(wf, rg, parent="commonLogicle")
xyplot(`SSC-A` ~ `FSC-A` , data = wf[["commonLogicle"]], subset = infection == "Cre",smooth=F,filter=rg,xlim=c(2,4.5),ylim=c(2,4.5))

# Ellipsoid gates

# A function to help covariance matrix calcculation for ellipsoid definition:

cov.matrix <- function (a, b, angle) {
   theta <- angle * (pi/180)
   c1 <- ((cos(theta)^2)/a^2) + ((sin(theta)^2)/b^2)
   c2 <- sin(theta) * cos(theta) * ((1/a^2) - (1/b^2))
   c3 <- ((sin(theta)^2)/a^2) + ((cos(theta)^2)/b^2)
   m1 <- matrix(c(c1, c2, c2, c3), byrow=TRUE, ncol=2)
   m2 <- solve(m1)
   m2
}

ellipsis<-c(130000,60000,35,100000,30000)

cov <- cov.matrix(ellipsis[4],ellipsis[5],ellipsis[3])
dimnames(cov)<-list(c("FSC-A", "SSC-A"), c("FSC-A", "SSC-A"))
eg <- ellipsoidGate(filterId= "LiveCells", .gate=cov, mean=c("FSC-A"=ellipsis[1], "SSC-A"=ellipsis[2]))
add(wf, eg, parent="Compensated Samples")

xyplot(`SSC-A` ~ `FSC-A` , data = wf[["Compensated Samples"]], subset = infection == "Cre",smooth=F,filter=eg)

#plot(Data(wf[["LiveCells+"]])[[1]],c("FSC-A","SSC-A"))

library(shape)
plot(Data(wf[["Compensated Samples"]])[[2]],c("FSC-A","SSC-A"),main="Compensated samples")
plotellipse(rx = ellipsis[4], ry = ellipsis[5], mid = ellipsis[1:2],angle=ellipsis[3],lcol="red")

par(mfrow=c(6,3), oma=c(0,0,0,0), mar=c(2,1,1,1))
for (i in 1:18){
  plot(Data(wf[["Compensated Samples"]])[[i]],c("FSC-A","SSC-A"),main=paste(pData(my.dat)$name[i],i,sep=" "))
  plotellipse(rx = ellipsis[4], ry = ellipsis[5], mid = ellipsis[1:2],angle=ellipsis[3],lcol="red")
}
par(mfrow=c(1,1))

#undo(wf)


ellipsis<-c(4.2,3.8,45,0.65,0.2)

cov <- cov.matrix(ellipsis[4],ellipsis[5],ellipsis[3])
dimnames(cov)<-list(c("FSC-A", "SSC-A"), c("FSC-A", "SSC-A"))
eg1 <- ellipsoidGate(filterId= "EllipsoidGate1", .gate=cov, mean=c("FSC-A"=ellipsis[1], "SSC-A"=ellipsis[2]))
add(wf, eg1, parent="commonLogicle")

xyplot(`SSC-A` ~ `FSC-A` , data = wf[["commonLogicle"]], subset = infection == "Cre",smooth=F,filter=eg1,xlim=c(2,4.5),ylim=c(2,4.5))


#plot(Data(wf[["EllipsoidGate1+"]])[[1]],c("FSC-A","SSC-A"))

par(mfrow=c(6,3), oma=c(0,0,0,0), mar=c(2,1,1,1))
for (i in 1:18){
  plot(Data(wf[["commonLogicle"]])[[i]],c("FSC-A","SSC-A"),main=paste(pData(my.dat)$name[i],i,sep=" "),xlim=c(2.9,4.6),ylim=c(2.5,4.6))
  plotellipse(rx = ellipsis[4], ry = ellipsis[5], mid = ellipsis[1:2],angle=ellipsis[3],lcol="red")
}
par(mfrow=c(1,1))


############################################################
# Polygon gate of living cells
############################################################

#undo(wf)

poly<-matrix(c(25000,0,220000,230000,300000,230000,300000,50000,220000,50000,100000,0),ncol=2,nrow=6,byrow=T)
colnames(poly) <- c("FSC-A","SSC-A")
pg <- polygonGate(filterId="LivingCellsPoly", .gate=poly)
add(wf, pg, parent="Compensated Samples")

par(mfrow=c(6,3), oma=c(0,0,0,0), mar=c(2,1,1,1))
for (i in 1:18){
  plot(Data(wf[["Compensated Samples"]])[[i]],c("FSC-A","SSC-A"),main=paste(pData(my.dat)$name[i],i,sep=" "))
  polygon(poly,border="red")
}
par(mfrow=c(1,1))

xyplot(`SSC-A` ~ `FSC-A` , data = wf[["Compensated Samples"]], subset = infection == "Cre",smooth=F,filter=pg)

#undo(wf)

poly<-matrix(c(3.7,3.6, 4.5,4.3, 4.6,4.3, 4.6,3.9, 4.5,3.7, 3.7,2.9),ncol=2,nrow=6,byrow=T)
colnames(poly) <- c("FSC-A","SSC-A")
pg2 <- polygonGate(filterId="LivingCellsPolyTransformed", .gate=poly)
add(wf, pg2, parent="commonLogicle")

par(mfrow=c(6,3), oma=c(0,0,0,0), mar=c(2,1,1,1))
for (i in 1:18){
  plot(Data(wf[["commonLogicle"]])[[i]],c("FSC-A","SSC-A"),main=paste(pData(my.dat)$name[i],i,sep=" "),xlim=c(2.8,4.5),ylim=c(2,4.6))
  polygon(poly,border="red")
}
par(mfrow=c(1,1))

xyplot(`SSC-A` ~ `FSC-A` , data = wf[["commonLogicle"]], subset = infection == "Cre",smooth=F,filter=pg2,xlim=c(2,4.7),ylim=c(2,4.5))

############################################################
# Gate CD4+ GFP+ cells
############################################################


rg2<-rectangleGate(filterId="RectangleGFP","GFP-A"=c(2.1, 3.5), "PerCP-Cy5-5-A"=c(2.5, 4))
add(wf, rg2, parent="EllipsoidGate1+")

par(mfrow=c(6,3), oma=c(0,0,0,0), mar=c(2,2,3,1))
for (i in 1:18){
  plot(Data(wf[["EllipsoidGate1+"]])[[i]],c("PerCP-Cy5-5-A","GFP-A"),main=paste(pData(my.dat)$name[i],i,sep=" "))
  mtext("CD4",side=1,line=1,cex=0.6)
  mtext("GFP",side=2,line=1,cex=0.6)
  rect(2.5,2.1,4,3.5,border="red")
}
par(mfrow=c(1,1),oma=c(1,1,1,1))


rg3<-rectangleGate(filterId="RectangleGFP1","GFP-A"=c(1000, 8000), "PerCP-Cy5-5-A"=c(3000, 22000))
add(wf, rg3, parent="LiveCells+")

par(mfrow=c(6,3), oma=c(0,0,0,0), mar=c(2,2,3,1))
for (i in 1:18){
  plot(Data(wf[["LiveCells+"]])[[i]],c("PerCP-Cy5-5-A","GFP-A"),main=paste(pData(my.dat)$name[i],i,sep=" "),xlim=c(0,30000),ylim=c(0,10000))
  mtext("CD4",side=1,line=1,cex=0.6)
  mtext("GFP",side=2,line=1,cex=0.6)
  rect(3000,1000,22000,8000,border="red")
}
par(mfrow=c(1,1),oma=c(1,1,1,1))


#undo(wf)
rg4<-rectangleGate(filterId="GFP+CD4+","GFP-A"=c(2.1, 3.5), "PerCP-Cy5-5-A"=c(2.5, 4))
add(wf, rg4, parent="LivingCellsPolyTransformed+")

par(mfrow=c(6,3), oma=c(0,0,0,0), mar=c(2,2,3,1))
for (i in 1:18){
  plot(Data(wf[["LivingCellsPolyTransformed+"]])[[i]],c("GFP-A","PerCP-Cy5-5-A"),main=paste(pData(my.dat)$name[i],i,sep=" "))
  mtext("GFP",side=1,line=1,cex=0.6)
  mtext("CD4",side=2,line=1,cex=0.6)
  rect(2.1,2.5, 3.5,4 ,border="red")
}
par(mfrow=c(1,1),oma=c(1,1,1,1))
xyplot(`PerCP-Cy5-5-A` ~ `GFP-A` , data = wf[["LivingCellsPolyTransformed+"]], subset = infection == "Cre",smooth=F,filter=rg4)


############################################################
# Use a quad gate for that
############################################################

qg<-quadGate(filterId="QuadCD_GFP", "GFP-A"=2.1, "PerCP-Cy5-5-A"=2.5)
add(wf, qg, parent="LivingCellsPolyTransformed+")
xyplot(`PerCP-Cy5-5-A` ~ `GFP-A` , data = wf[["LivingCellsPolyTransformed+"]], subset = infection == "Cre",smooth=F,filter=qg)


############################################################
# Gate IL2+ cells
############################################################

#undo(wf)
rg5<-rectangleGate(filterId="IL2","PE-A"=c(1.7, 4), "PerCP-Cy5-5-A"=c(2.4, 4.1))
add(wf, rg5, parent="GFP+CD4+")

par(mfrow=c(6,3), oma=c(0,0,0,0), mar=c(2,2,3,1))
for (i in 1:18){
  plot(Data(wf[["RectangleGFP0+"]])[[i]],c("PE-A","PerCP-Cy5-5-A"),main=paste(pData(my.dat)$name[i],i,sep=" "),ylim=c(1.6,4.2))
#  contour(Data(wf[["RectangleGFP0+"]])[[i]],c("PE-A","PerCP-Cy5-5-A"),main=paste(pData(my.dat)$name[i],i,sep=" "),ylim=c(1.6,4.2))
  mtext("IL2",side=1,line=1,cex=0.6)
  mtext("CD4",side=2,line=1,cex=0.6)
  rect(1.7,2.4, 4,4.1 ,border="red")
}
par(mfrow=c(1,1),oma=c(1,1,1,1))

xyplot(`PerCP-Cy5-5-A` ~ `PE-A` , data = wf[["GFP+CD4++"]], subset = infection == "Cre",smooth=F,filter=rg5)
xyplot(`PerCP-Cy5-5-A` ~ `PE-A` , data = wf[["GFP+CD4++"]], subset = infection == "Migr",smooth=F,filter=rg5)

# plot histograms

par(mfrow=c(6,3), oma=c(0,0,0,0), mar=c(2,2,3,1))
for (i in 1:18){
  plot(Data(wf[["IL2+"]])[[i]],"PE-A",main=paste(pData(my.dat)$name[i],i,sep=" "),ylim=c(0,300))
  mtext("IL2",side=2,line=1,cex=0.6)
}
par(mfrow=c(1,1),oma=c(1,1,1,1))

############################################################
# Automatic filtering of IL2+ cells
############################################################

ag<-kmeansFilter("PE-A"=c("kMeanIL2-","kMeanIL2+"), filterId="Il2KmFilter")
add(wf, ag, parent="GFP+CD4+")
summary(wf[["IL2+"]])
summary(wf[["kMeanIL2+"]])

############################################################
# Statistics
############################################################

il2results<-cbind(pData(my.dat),summary(wf[["IL2+"]]))

cre.res<-data.frame(Stim.00h=il2results[il2results$infection=="Cre" & il2results$stimulation=="00h","percent"],Stim.05h=il2results[il2results$infection=="Cre" & il2results$stimulation=="05h","percent"],Stim.24h=il2results[il2results$infection=="Cre" & il2results$stimulation=="24h","percent"])

migr.res<-data.frame(Stim.00h=il2results[il2results$infection=="Migr" & il2results$stimulation=="00h","percent"],Stim.05h=il2results[il2results$infection=="Migr" & il2results$stimulation=="05h","percent"],Stim.24h=il2results[il2results$infection=="Migr" & il2results$stimulation=="24h","percent"])

barplot(as.matrix(cre.res),beside=T,main="Cre infection",ylab="%")
barplot(as.matrix(migr.res),beside=T,main="Migr infection",ylab="%")

boxplot(percent ~ stimulation * infection,data=il2results,ylab="%",main="IL2 expressing cells",names=c("Cre 0h","Cre 5h","Cre 24h","Migr 0h","Migr 5h","Migr 24h"),col=c(rep("red",3),rep("blue",3)))


boxplot(percent ~ infection * stimulation,data=il2results,ylab="%",main="IL2 expressing cells",names=c("Cre 0h","Migr 0h","Cre 5h","Migr 5h","Cre 24h","Migr 24h"),col=c("red","blue"))

############################################################
# The same for the automated gate
############################################################

il2results.ag<-cbind(pData(my.dat),summary(wf[["kMeanIL2+"]]))

cre.res.ag<-data.frame(Stim.00h=il2results.ag[il2results.ag$infection=="Cre" & il2results.ag$stimulation=="00h","percent"],Stim.05h=il2results.ag[il2results.ag$infection=="Cre" & il2results.ag$stimulation=="05h","percent"],Stim.24h=il2results.ag[il2results.ag$infection=="Cre" & il2results.ag$stimulation=="24h","percent"])

migr.res.ag<-data.frame(Stim.00h=il2results.ag[il2results.ag$infection=="Migr" & il2results.ag$stimulation=="00h","percent"],Stim.05h=il2results.ag[il2results.ag$infection=="Migr" & il2results.ag$stimulation=="05h","percent"],Stim.24h=il2results.ag[il2results.ag$infection=="Migr" & il2results.ag$stimulation=="24h","percent"])

barplot(as.matrix(cre.res.ag),beside=T,main="Cre infection",sub="automated gating",ylab="%")
barplot(as.matrix(migr.res.ag),beside=T,main="Migr infection",sub="automated gating",ylab="%")

boxplot(percent ~ stimulation * infection,data=il2results.ag,ylab="%",main="IL2 expressing cells",names=c("Cre 0h","Cre 5h","Cre 24h","Migr 0h","Migr 5h","Migr 24h"),col=c(rep("red",3),rep("blue",3)),sub="automated gating")


boxplot(percent ~ infection * stimulation,data=il2results.ag,ylab="%",main="IL2 expressing cells",names=c("Cre 0h","Migr 0h","Cre 5h","Migr 5h","Cre 24h","Migr 24h"),col=c("red","blue"))



t.test(migr.res$Stim.00h,migr.res$Stim.05h,alternative="greater")














############################################################
# Alternatives
############################################################


# More advanced transformations with flowTrans

library(flowTrans)

#tr1<-lapply(as(my.dat,"list"),function(x)flowTrans(dat=x,fun="mclMultivArcSinh",dims=c("FSC-A","SSC-A","GFP-A","PE-A","PerCP-Cy5-5-A"),n2f=FALSE,parameters.only=FALSE));

#transformed2<-flowTrans(dat=my.dat,fun="mclMultivArcSinh",dims=c("FSC-H","SSC-H"),n2f=FALSE,parameters.only=TRUE)


#add(wf, tr1, parent = "Compensated Samples")

#plot(my.dat[[1]], c("FSC-A", "SSC-A"),main="Untransformed")
#X11()
#plot(tr1[[1]]$result, c("FSC-A", "SSC-A"),main="ArcSinh transformation")


# Box-cox for entire Set

my.dat.boxcox<-lapply(as(my.dat,"list"),function(x)flowTrans(dat=x,fun="mclMultivBoxCox",dims=c("FSC-A","SSC-A","GFP-A","PE-A","PerCP-Cy5-5-A"),n2f=FALSE,parameters.only=FALSE));



#add(wf, my.dat.boxcox, parent = "Compensated Samples")


par(mfrow=c(3,2))

for (i in 1:3){
  plot(my.dat[[i]], c("FSC-A", "SSC-A"),sub="Untransformed",main=names(tr2)[i])
  plot(my.dat.boxcox[[i]]$result, c("FSC-A", "SSC-A"),sub=paste("Box-Cox transformation Theta: ",extractParams(my.dat.boxcox[[i]])$'FSC-A'),main=names(tr2)[i],ylim=c(10,50))

}


#for (i in 1:3){
#  plot(my.dat[[i]],main="Untransformed")
#  X11()
#  plot(my.dat.boxcox[[i]]$result,main="BoxCox transformation",sub=paste("Theta: ",extractParams(my.dat.boxcox[[i]])$'FSC-A'))
#  X11()
#}
