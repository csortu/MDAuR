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
