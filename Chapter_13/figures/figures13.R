### Figures for chapter 13
##library(Cairo)

## figure 13.1

GFP.exc<-read.csv("EGFPpH7_exc.csv")
GFP.em<-read.csv("EGFPpH7_em.csv")

tiff("figure13_1.tif",width=5,height=5,units="in",res=1200)
par(cex=0.8)

plot(GFP.exc,
     xlim=c(450,550),type="l",
     main="Excitation and emission spectrum of GFP",
     xlab="Wavelength [nm]",
     ylab="",ylim=c(0,110),yaxt="n")
points(GFP.em,type="l",lty=5)
abline(v=488,lwd=2,lty=3)

legend(445,115,c("Excitation","Emission","Argon 488 Laser"),lty=c(1,5,3),bty="n")

arrows(490,105,508,105,code=3,length=0.15)
text(500,110,"Wavelength shift",adj=0.5,cex=0.8)

dev.off()

## figure 13.2

GFP.exc<-read.csv("EGFPpH7_exc.csv")
GFP.em<-read.csv("EGFPpH7_em.csv")

PE.exc<-read.csv("R-PE-801ph75_exc.csv")
PE.em<-read.csv("R-PE-801ph75_em.csv")

tiff("figure13_2.tif",width=5,height=5,units="in",res=1200)
par(cex=0.8,mfrow =c(2,1),mai=c(0.52,0.67,0.42,0.12))

plot(NA, 
     ylim=c(0,110),xlim=c(250,650), 
     type="l", 
     col="red",
     main="GFP", xlab="Wavelength [nm]", ylab="",
     yaxt="n")
#polygon(c(495,495,525,525),c(0,100,100,0),col="gray")
rect(495,0,525,100,density = 20)
polygon(seq(545,625),c(0,GFP.em$em[GFP.em$wl>545 & GFP.em$wl< 625],0),col="gray")

points(GFP.exc,type="l")
points(GFP.em,type="l",lty=5)
abline(v=488,lwd=2,lty=3)
legend(310,115,c("Excitation","Emission","Argon 488 Laser","Detector"),lty=c(1,5,3,0),bty="n",cex=0.8)
rect(322,62,347,72,density = 20)
mtext("A", side = 3, line = 0, adj = 0, cex = 1.3)
text(585,60,labels = "spilling",cex=0.7,adj=0.5)
arrows(585,55,575,20,length = 0.05)

plot(NA,
     ylim=c(0,110),xlim=c(250,650),
     type="l",col="red",
     main="R-PE",
     xlab="Wavelength [nm]",ylab="",
     yaxt="n")

#polygon(c(545,545,625,625),c(0,100,100,0),col="gray")
rect(545,0,625,100,density = 20)

points(PE.exc,type="l")
points(PE.em,type="l",lty=5)
abline(v=488,lty=3,lwd=2)
mtext("B", side = 3, line = 0, adj = 0, cex = 1.3)

dev.off()

## figure 13.3

library(flowCore)
library(flowViz)
library(plotrix)

fcs<-read.FCS("Cells_300 000.fcs")

my.dat <- exprs(fcs)[,1:2]
my.dat <- my.dat/10

tiff("figure13_3.tif",width=5,height=5,units="in",res=1200)

plot(my.dat[1:100000,],
     pch=20,
     cex=0.2,main="Human peripherial blood cells",
     col=gray(0.2))
draw.ellipse(9000,14000,6000,3000,angle = 100,lwd=3)
draw.ellipse(17000,6000,4000,2000,angle = 0,lwd=3)
draw.ellipse(13000,2700,4000,1000,angle = 15,lwd=3)

text(12000,22000,labels = "Granulocytes")
text(18000,10000,labels = "Monocytes")
text(9000,5500,labels = "Lymophocytes")

dev.off()

## figure 13.4

library(flowCore)
library(flowViz)
fcspath<-"" # set the path where the fcs_files.zip is opened
fcs.files<-list.files(fcspath,pattern=".fcs",full=T)
dat<-read.flowSet(fcs.files)
my.dat<-dat[5:22]
max.fsc<-max(exprs(my.dat[[1]])[,"FSC-A"])
max.ssc<-max(exprs(my.dat[[1]])[,"SSC-A"])
margin.cells<-as.logical((exprs(my.dat[[1]])[,"SSC-A"]==max.ssc) + (exprs(my.dat[[1]])[,"FSC-A"]==max.fsc))

tiff("figure13_4.tif",width=5,height=5,units="in",res=1200)
plot(my.dat[[1]],c("FSC-A","SSC-A"),smooth=F,main="Margin events")
points(exprs(my.dat[[1]])[margin.cells,c("FSC-A","SSC-A")])
dev.off()

## compensation in workflow
library(flowCore)
library(flowViz)
fcspath<-"" # set the path where the fcs_files.zip is opened
fcs.files<-list.files(fcspath,pattern=".fcs",full=T)
dat<-read.flowSet(fcs.files)
pData(dat)<-data.frame(name=pData(dat)[, 1],infection=NA,stimulation=NA,row.names=row.names(pData(dat)),desc=NA)
varMetadata(phenoData(dat))$labelDescription <- c("File","Infection","Stimulation","Description")
pData(dat)$infection<-c(rep(NA,4),rep("Cre",9),rep("Migr",9),NA)
pData(dat)$stimulation<- c(rep(NA,4),rep(c("24h","05h","00h"),6),NA)
pData(dat)$desc<-fsApply(dat,function(f) f@description $`TUBE NAME`)
my.dat<-dat[5:22]

for (i in 1:18){
  max.fsc<-max(exprs(my.dat[[i]])[,"FSC-A"])
  max.ssc<-max(exprs(my.dat[[i]])[,"SSC-A"])
  margin.cells<-as.logical((exprs(my.dat[[i]])[,"SSC-A"]==max.ssc) + (exprs( my.dat[[i]])[,"FSC-A"]==max.fsc))
  print(pData(my.dat)$name[i])
  print(paste("Margin cell ratio",sum(margin.cells)/nrow(my.dat[[i]])))
  my.dat[[i]]<-(my.dat[[i]])[!margin.cells]
}


wf <- workFlow(my.dat)
add(wf, compensation(spill, compensationId = "Compensated Samples"))

tiff("figure13_5_color.tif",width=5,height=5,units="in",res=1200)
xyplot(`SSC-A` ~ `FSC-A` , 
       data = wf[["Compensated Samples"]], 
       subset = infection == "Cre",
       main="Compensated samples: Cre infection",
       strip=strip.custom( par.strip.text=list(cex=.5)))
dev.off()

## gating

## ellypsis gate

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

tiff("figure13_6_color.tif",width=5,height=5,units="in",res=1200)

xyplot(`SSC-A` ~ `FSC-A` , 
       data = wf[["Compensated Samples"]], 
       subset = infection == "Cre",
       smooth=F,filter=eg,
       main="Living cells in Cre samples with manual gating",
       strip=strip.custom( par.strip.text=list(cex=.5)))

dev.off()

## poly gate

poly<-matrix(c(25000,0,220000,230000,300000,230000,300000,50000,220000,50000,100000,0),ncol=2,nrow=6,byrow=T)
colnames(poly) <- c("FSC-A","SSC-A")
pg <- polygonGate(filterId="LivingCellsPoly", .gate=poly)
add(wf, pg, parent="Compensated Samples")

tiff("figure13_7_color.tif",width=5,height=5,units="in",res=1200)

xyplot(`SSC-A` ~ `FSC-A` , 
       data = wf[["Compensated Samples"]], 
       subset = infection == "Cre",
       smooth=F,filter=pg,
       main="Living cells in Cre samples with polygon gating",
       strip=strip.custom( par.strip.text=list(cex=.5)))

dev.off()

## statistics

rm(wf)
wf <- workFlow(my.dat)
add(wf, compensation(spill, compensationId = "Compensated Samples"))
tr<-estimateLogicle(flowFrame(fsApply(my.dat, function(x) exprs(x))),colnames(my.dat)[1:5])
identifier(tr) <- "commonLogicle"
add(wf, tr, parent = "Compensated Samples")

poly<-matrix(c(3.7,3.6, 4.5,4.3, 4.6,4.3, 4.6,3.9, 4.5,3.7, 3.7,2.9),ncol=2,nrow=6,byrow=T)
colnames(poly) <- c("FSC-A","SSC-A")
pg2 <- polygonGate(filterId="LivingCellsPolyTransformed", .gate=poly)
add(wf, pg2, parent="commonLogicle")

rg4<-rectangleGate(filterId="GFP+CD4+","GFP-A"=c(2.1, 3.5), "PerCP-Cy5-5-A"=c(2.5, 4))
add(wf, rg4, parent="LivingCellsPolyTransformed+")

rg5<-rectangleGate(filterId="IL2","PE-A"=c(1.7, 4), "PerCP-Cy5-5-A"=c(2.4, 4.1))
add(wf, rg5, parent="GFP+CD4++")

il2results<-cbind(pData(my.dat),summary(wf[["IL2+"]]))

tiff("figure13_8.tif",width=5,height=5,units="in",res=1200)

boxplot(percent ~ infection * stimulation,
        data=il2results,
        ylab="%",main="IL2 expressing cells",
        names=c("Cre 0h","Migr 0h","Cre 5h","Migr 5h","Cre 24h","Migr 24h"),
        col=c(gray(1:2/3)),cex.axis=0.7)

dev.off()
