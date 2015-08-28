GFP.exc<-read.csv("EGFPpH7_exc.csv")
GFP.em<-read.csv("EGFPpH7_em.csv")

plot(GFP.exc,xlim=c(250,650),type="l",col="red",main="Fluorescent excitation and emission spectrum of GFP",xlab="Wavelength [nm]",ylab="")
points(GFP.em,type="l",col="blue")
abline(v=488,col="yellow",lwd=2)

legend(300,90,c("Excitation","Emission","Argon 488 Laser"),col=c("red","blue","yellow"),lty=1)

# excitation peak at 490 nm 
# emission peak at 509 nm

# Common laser to excite: 488 nm
# Common filter to measure: 525 nm +- 25

PE.exc<-read.csv("R-PE-801ph75_exc.csv")
PE.em<-read.csv("R-PE-801ph75_em.csv")

plot(NA,ylim=c(0,100),xlim=c(250,650),type="l",col="red",main="Fluorescent excitation and emission spectrum of GFP an R-PE",xlab="Wavelength [nm]",ylab="")

polygon(c(495,495,525,525),c(0,100,100,0),col="green")
polygon(c(545,545,625,625),c(0,100,100,0),col="orange")
polygon(seq(545,625),c(0,GFP.em$em[GFP.em$wl>545 & GFP.em$wl< 625],0),col="gray")

points(GFP.exc,type="l",col="red")
points(GFP.em,type="l",col="blue")
points(PE.exc,type="l",col="red",lty=2)
points(PE.em,type="l",col="blue",lty=2)
abline(v=488,col="red",lty=4,lwd=2)

legend(300,90,c("GFP Excitation","GFP Emission","PE Excitation","PE Emission","Argon 488 Laser"),col=c("red","blue","red","blue","red"),lty=c(1,1,2,2,4),bty="n")
legend(310,70,c("GFP filter","PE filter","GFP spill \nto PE channel"),fill=c("green","orange","gray"),bty="n")
