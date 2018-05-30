## Required Packages
library("betapart")
library("reshape")
library("vegan")
library("flowWorkspace")
library("ggplot2")
library("ggcyto")
library("plyr")
library("car")
library("corrplot")
library("pangaear")
library("ggbiplot")
library("psych")
library("cowplot")
library("RColorBrewer")
library("dunn.test")
library("seriation")
library("ggbiplot")
library("gvlma")
library("RVAideMemoire")
library("gplots")
library("flowDiv")
## Workings directly
setwd("/home/brunomsw/Documentos/DOUTORADO/DOC2018/MAIN/100LAGOS/FCS")
ids=read.csv("../METADATA/ID_DILUTION.csv", sep="\t")
metadata=read.csv("../METADATA/MF.csv", dec=",")

cem_lagos = flowDiv("Bruno_100Lagos_Beadsv10.wsp", "Bact")

# Limpeza das planilhas (remoção de dados excedentes)
lakes<-unlist(lapply(names(cem_lagos$Alpha), function(x)substr(x, 1, 23)))
ids=ids[ids$FCS%in%lakes,]
metadata=metadata[metadata$lakeID%in%ids$LAGO,]
#
match(metadata$lakeID, ids$LAGO)
metadata=metadata[match(ids$LAGO, metadata$lakeID),]
metadata=metadata[,complete.cases(t(metadata))]

# Update
cem_lagos = flowDiv("Bruno_100Lagos_Beadsv10.wsp", "Bact", dilutions = ids$DILUTION)

# mat1=as.dist(cem_lagos$Beta) Nao funciona bem com envfit...
mat2<- metaMDS(cem_lagos$Matrices, autotransform = F, try=200)
envfit(mat2, metadata[,c(8:73, 84:87)])
#


b=gvlma(lm(cem_lagos$Alpha~metadata$PB.µgC.l.h.)) #OK

# Função para verificar a adequação (gvlma OK pra todos os testes)
gvbruno<-function(a){
soma=sum(a$GlobalTest$GlobalStat4$Decision, a$GlobalTest$DirectionalStat1$Decision, a$GlobalTest$DirectionalStat2$Decision, a$GlobalTest$DirectionalStat3$Decision, a$GlobalTest$DirectionalStat4$Decision)
return(soma==0)}

# Testando
for(i in 8:73) print(gvbruno(gvlma(lm(cem_lagos$Alpha~metadata[,i]))))
for(i in 8:73) print(gvbruno(gvlma(lm(cem_lagos$Pielou~metadata[,i]))))

# Correlaçao
a=c(); b=c()
for(i in 8:73) a[i-7]=cor.test(cem_lagos$Alpha, metadata[,i], method = "spearman", exact = F)$p.value
for(i in 8:73) b[i-7]=cor.test(cem_lagos$Pielou, metadata[,i], method = "spearman", exact = F)$p.value
cbind(a,b)


snoop<-cbind(a,b, colnames(metadata)[8:73])
# Alpha
unlist(snoop[,3][snoop[,1]<=0.05])
#Pielou
unlist(snoop[,3][snoop[,2]<=0.05])


# Matriz renomeada
m1=cem_lagos$Matrices

#Definindo os indices
rich=apply(m1, 1, function(x) specnumber(x, MARGIN = 1))
pielou=apply(m1, 1, function(x)diversity(x)/log(vegan::specnumber(x)))
shannon=apply(m1, 1, function(x) diversity(x,index = "shannon"))
simpson=apply(m1, 1, function(x) diversity(x,index = "simpson"))
invsimpson=apply(m1, 1, function(x) diversity(x,index = "invsimpson"))

dataall<-cbind(rich, pielou, shannon, simpson, invsimpson,metadata[,8:73])
par(mfrow=c(4,4))
for (i in 1:71) plot(dataall$rich, dataall[,i], ylab=names[i])
for (i in 1:71) plot(dataall$shannon, dataall[,i], ylab=names[i])
for (i in 1:71) plot(dataall$pielou, dataall[,i], ylab=names[i])


# Correlaçao 2
a=c(); b=c(); d=c(); e=c(); f=c()
for(i in 8:73) a[i-7]=cor.test(rich, metadata[,i], method = "spearman", exact = F)$p.value
for(i in 8:73) b[i-7]=cor.test(pielou, metadata[,i], method = "spearman", exact = F)$p.value
for(i in 8:73) d[i-7]=cor.test(shannon, metadata[,i], method = "spearman", exact = F)$p.value
for(i in 8:73) e[i-7]=cor.test(simpson, metadata[,i], method = "spearman", exact = F)$p.value
for(i in 8:73) f[i-7]=cor.test(invsimpson, metadata[,i], method = "spearman", exact = F)$p.value


cbind(a,b,d, e, f)


snoop<-cbind(a,b,d, e,f, colnames(metadata)[8:73])
names=c("a","d","d","e","f",colnames(metadata)[8:73])
# Rich
unlist(snoop[,6][snoop[,3]<=0.05])
# Shannon
unlist(snoop[,6][snoop[,4]<=0.05])
# Simp
unlist(snoop[,6][snoop[,5]<=0.05])
# Inv
unlist(snoop[,6][snoop[,6]<=0.05])



