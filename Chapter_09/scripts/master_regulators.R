library(Biobase)
load("eset_tf.rda")

eset.tf
plot(density(log2(exprs(eset.tf)),na.rm=T))
featureNames(eset.tf)

# probes for key t-cell TFs
# 1341_at - SPI1 = PU.1
# 40511_at - GATA3
# 41504_s_at 41505_r_at - MAF
# TBX21 = T-bet is not on the chip!
# 33592_at - ZBTB7B = Th-POK


target.tf.probes<-c("1341_at","40511_at","41504_s_at","33592_at")
#target.tf.probes<-c("1341_at","40511_at","41504_r_at","33592_at")
names(target.tf.probes)<-c("SPI1","GATA3","MAF","ZBTB7B")


#library(hgu95av2.db)
#annot<-select(hgu95av2.db,keys=featureNames(eset.tf),keytype="PROBEID",columns=c("PROBEID","ENTREZID","SYMBOL"))
#annot$ENTREZID<-as.numeric(annot$ENTREZID)
#write.csv(annot,"annotation.csv")

library(RTN)

annot2<-read.csv("annot.csv",row.names=1)
tf.rtni<-new("TNI", gexp=exprs(eset.tf),transcriptionFactors=target.tf.probes)
tf.rtni.pp<-tni.preprocess(tf.rtni,gexpIDs=annot2)
tf.rtni.per<-tni.permutation(tf.rtni.pp,estimator='kendall',pValueCutoff=0.01)
tf.rtni.boot<-tni.bootstrap(tf.rtni.per,estimator='kendall',consensus=95)
tf.rtni.dpi<-tni.dpi.filter(tf.rtni.boot)

tni.get(tf.rtni.dpi,what="summary")


tf.rtni<-new("TNI", gexp=exprs(eset.tf),transcriptionFactors=target.tf.probes)
tf.rtni<-tni.preprocess(tf.rtni,gexpIDs=annot2)
tf.rtni<-tni.permutation(tf.rtni,estimator='kendall',pValueCutoff=0.01)
tf.rtni<-tni.bootstrap(tf.rtni,estimator='kendall',consensus=95)
tf.rtni<-tni.dpi.filter(tf.rtni)

tni.get(tf.rtni,what="summary")

g<-tni.graph(tf.rtni)
V(g)$label<-as.character(annot2[annot2$PROBEID %in% V(g)$name,"SYMBOL"])
plot(g,vertex.size=20)

