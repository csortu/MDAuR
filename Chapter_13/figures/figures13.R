### Figures for chapter 13
##library(Cairo)

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


