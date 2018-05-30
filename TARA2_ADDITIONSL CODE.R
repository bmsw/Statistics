mOTU_abs=read.csv("/home/brunomsw/Documentos/DOUTORADO/DATA_FCS/Tara/TARA METADATA/miTAG.taxonomic.profiles.release.tsv.csv", sep="\t", check.names = F)
mOTU=read.csv("/home/brunomsw/Documentos/DOUTORADO/DATA_FCS/Tara/TARA METADATA/mOTU.linkage-groups.relab.release.csv", sep="\t", check.names = F)
eggNOG_Proc = read.csv("/home/brunomsw/Documentos/DOUTORADO/DATA_FCS/Tara/TARA METADATA/TARA139.og.prok.fpkm.release.csv", sep="\t", check.names = F)
eggNOG=read.csv("/home/brunomsw/Documentos/DOUTORADO/DATA_FCS/Tara/TARA METADATA/TARA243.OG.profile.release.csv", sep="\t", check.names = F)

eggNOG=eggNOG[-1,]
eggNOG_Proc = eggNOG_Proc[-63771,]
mOTU=mOTU[-1,]
mOTU_abs=mOTU_abs[-1,]

colnames(eggNOG)[80]="TARA_148_MES_0.22-3" 
colnames(eggNOG_Proc)[52]="TARA_148_MES_0.22-3"
colnames(mOTU)[88]="TARA_148_MES_0.22-3"
colnames(mOTU_abs)[94]="TARA_148_MES_0.22-3"






newds <-function(x){
  
  t1=x[,c(grep("0.22-3",colnames(x)))]
  t2= t1[,-grep("MIX",colnames(t1), ignore.case = T)]
  return(t2)
}




mantel.fd<-function(x, y){
  
  d1=newds(x)
  m2=m1
  rownames(m2) = paste(rownames(m2), "_0.22-3", sep="")
  match1=match(rownames(m2), colnames(x))
  
  
  d2=x[,match1]
  d4=vegdist(t(d2), y)
  d3=vegdist(m1,method = "jaccard", binary = T)
 # d3=vegdist(m1)
  
  d31=  upperTriangle(as.matrix(d3))
  d41=  upperTriangle(as.matrix(d4))
  
  d5=data.frame(cbind(d31, d41))
  g1=ggplot(aes(d31, d41), data=d5)+geom_point(size=.5)+ geom_smooth(method='lm', colour="red")
  
  #plot(d31,d41, pch=".")
  #abline(lm(d41 ~ d31), col="red")
  print(cor(d31, d41))
  
  #return(mantel(d3, d4))
  return(g1)
}




mantel.fd(eggNOG, "bray")
mantel.fd(eggNOG, "euclidian")
mantel.fd(eggNOG_Proc, "bray")
mantel.fd(eggNOG_Proc, "euclidian")
mantel.fd(mOTU, "bray")
mantel.fd(mOTU, "euclidian")
mantel.fd(mOTU_abs, "bray")
mantel.fd(mOTU_abs, "euclidian")  







