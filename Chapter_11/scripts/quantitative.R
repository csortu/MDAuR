library(MSnbase)

mztab <- 'F063721.dat-mztab.txt'
qnt <- readMzTabData(mztab, what = "PEP")


sampleNames(qnt)<-c("TMT6.126","TMT6.127","TMT6.128","TMT6.129","TMT6.130","TMT6.131")

qntS <- normalise(qnt, "sum")
qntV <- normalise(qntS, "vsn")
qntV2 <- normalise(qnt, "vsn")

par(mfrow=c(2,2))
boxplot(exprs(qnt),main="Raw data")
boxplot(exprs(qntS),main="Sum normalization")
boxplot(exprs(qntV),main="Variance stabilization of sum")
boxplot(exprs(qntV2),main="Variance stabilization on raw data")
par(mfrow=c(1,1))


protqnt <- combineFeatures(qnt, groupBy = fData(qnt)$accession, fun = sum)

protexp<-exprs(protqnt)

row.names(protexp)[400:404]<-c("PYGM_RABIT","TRYP_PIG","ENO1_YEAST","ALBU_BOVIN","CYC_BOVIN")

boxplot(protexp[1:399,])

matplot(t(protexp[400:404,]), type = "b",col=c("red","blue","black","green","cyan"),xlab="Samples",ylab="Protein abundancy")
legend("topright",row.names(protexp)[400:404],lty = 1, bty = "n", cex = .8, col = c("red","blue","black","green","cyan"))

wbcol <- colorRampPalette(c("white", "darkblue"))(256)
heatmap(protexp[350:404,],col=wbcol)

