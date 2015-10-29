library(hpar)

id<-"ENSG00000158813"

eda.tissue<-getHpa(id, "NormalTissue")
eda.Rna<-getHpa(id, "Rna")
eda.loc<-getHpa(id, "SubcellularLoc")

