## Required Packages

library("reshape")

library("ggplot2")
library("ggcyto")
library("plyr")
library("car")
library("corrplot")
library("pangaear")
library("ggbiplot")

library("cowplot")
library("RColorBrewer")
library("dunn.test")
library("seriation")
library("ggbiplot")
library("gvlma")
library("RVAideMemoire")
library("gplots")

#####################
library("vegan")
library("flowWorkspace")
library("betapart")
####################### Sourcing codes ##################################
source("/home/bruno/Documentos/DOC/flowDiv_SUBMISSION/Packages_and_Codes/flowDiv2/R/flowDiv_internalsCORRECTED.R")
source("/home/bruno/Documentos/DOC/flowDiv_SUBMISSION/Packages_and_Codes/flowDiv2/R/flowDiv_main2_2018teste.R")

## Diretorio de trabalho e Metadados
setwd("/home/bruno/Documentos/DOC/MAIN/100LAGOS/FCS")
ids=read.csv("../METADATA/ID_DILUTION.csv", sep="\t")
metadata=read.csv("../METADATA/MF.csv", dec=",")
metadata=metadata[,-c(26:57)]

# rimeira: registrar ID's para organização das linhas dos metadados
cem_lagos = flowDiv("Bruno_100Lagos_Beadsv10.wsp", gate.name = "Bact", beads = "Beads", use.beads = T )
2 3 5 
43
# Limpeza das planilhas: Alguns lagos não constam nas planilhas
lakes<-unlist(lapply(names(cem_lagos$Alpha), function(x)substr(x, 1, 23)))
ids=ids[ids$FCS%in%lakes,]
metadata=metadata[metadata$lakeID%in%ids$LAGO,]
# Reorganizando
match(metadata$lakeID, ids$LAGO)
metadata=metadata[match(ids$LAGO, metadata$lakeID),]
metadata=metadata[,complete.cases(t(metadata))]

# flowDiv com as diluições
cem_lagos = flowDiv("Bruno_100Lagos_Beadsv10.wsp", "Bact", dilutions = ids$DILUTION)
2 3 5
43
metadata=metadata[-94,] # Riacho da Cruz
m1=cem_lagos$Matrices[-94,-94] # Riacho da Cruz

library(FD)
traits<- data.frame(cem_lagos$ids[-94,1:3])
abund<-t(m1)

#
abund2=abund[!rowSums(abund)==0,]
traits2=traits[!rowSums(abund)==0,]

teste=dbFD(traits2, t(abund2))




# Definindo novos indices

dataall<-cbind(Ric=teste$FRic, Qrich=teste$qual.FRic, Ev=teste$FEve, Rao=teste$RaoQ, metadata[,c(4, 5, 8:41, 52:55)])

library("psych")

corr.test(dataall, method = "spearman")
res1 <- cor.mtest(dataall, 0.95)

library(corrplot)
corrplot(cor(dataall, method = "spearman", use="complete.obs"), p.mat = res1[[1]], sig.level=.05, type = "upper", pch.col="gray20", method = "square",
         col=c("red3", "blue3"), tl.col = "black", tl.cex = .01, pch.cex = 0.01)
################
a=cor(dataall, method = "spearman", use="complete.obs")[,1:4]
b=round(res1[[1]][,1:4], 2)
 cbind(a,b)

# Analise exploratoria
## Plots pareados
par(mfrow=c(3,3))
for (i in 6:46) plot(dataall$Ric, dataall[,i], ylab=colnames(dataall)[i])
for (i in 6:46) plot(dataall$Ev, dataall[,i], ylab=colnames(dataall)[i])


for (i in 6:47) plot(dataall$shannon, dataall[,i], ylab=names[i])
for (i in 6:47) plot(dataall$pielou, dataall[,i], ylab=names[i])

## Função para verificar a adequação (gvlma OK pra todos os testes)
gvbruno<-function(a){
  soma=sum(a$GlobalTest$GlobalStat4$Decision, a$GlobalTest$DirectionalStat1$Decision, a$GlobalTest$DirectionalStat2$Decision, a$GlobalTest$DirectionalStat3$Decision, a$GlobalTest$DirectionalStat4$Decision)
  return(soma==0)}

# Testando gvlma
gr=c(); gs=c(); gp=c(); gsi=c(); ginv=c()
for(i in 8:47) gr[i-7]=gvbruno(gvlma(lm(rich~metadata[,i])))
for(i in 8:47) gp[i-7]=gvbruno(gvlma(lm(pielou~metadata[,i])))
for(i in 8:47) gs[i-7]=gvbruno(gvlma(lm(shannon~metadata[,i])))
for(i in 8:47) gsi[i-7]=gvbruno(gvlma(lm(simpson~metadata[,i])))
for(i in 8:47) ginv[i-7]=gvbruno(gvlma(lm(invsimpson~metadata[,i])))


# Testando Correlações

r=c(); p=c(); sh=c(); si=c(); inv=c()
for(i in 8:47) r[i-7]=cor.test(rich, metadata[,i])$p.value
for(i in 8:47) p[i-7]=cor.test(pielou, metadata[,i])$p.value
for(i in 8:47) sh[i-7]=cor.test(shannon, metadata[,i])$p.value
for(i in 8:47) si[i-7]=cor.test(simpson, metadata[,i])$p.value
for(i in 8:47) inv[i-7]=cor.test(invsimpson, metadata[,i])$p.value

correls=cbind(cbind(p, r, sh, si, inv)<=0.05, gp, gr, gs,gsi, ginv, names[8:47])
correls<-data.frame(correls)

correls$p<-as.logical(correls$p)
correls$r<-as.logical(correls$r)
correls$sh<-as.logical(correls$sh)
correls$inv<-as.logical(correls$inv)
correls$si<-as.logical(correls$si)

correls$gp<-as.logical(correls$gp)
correls$gr<-as.logical(correls$gr)
correls$gs<-as.logical(correls$gs)
correls$gsi<-as.logical(correls$gsi)
correls$ginv<-as.logical(correls$ginv)


snoop1=cbind(p=(correls$p)&(correls$gp), r=(correls$r)&(correls$gr), sh=(correls$sh)&(correls$gs), si=(correls$si)&(correls$gsi), inv=(correls$inv)&(correls$ginv))
snoop2=cbind(snoop1, names[8:47])

snoop3=snoop2[!apply(snoop1,1, sum)==0,]

par(mfrow=c(3,3))
for (i in 1:19)plot(dataall$pielou, dataall[,snoop3[i,6]], ylab=snoop3[i,6])

# Texto
par(mfrow=c(3,3))
for (i in 1:19){plot(dataall$pielou, dataall[,snoop3[i,6]], ylab=snoop3[i,6], type="n"); text(dataall$pielou, dataall[,snoop3[i,6]], labels =metadata$lakeID)}











# Envfit traalha melhor com objetos metaMDS - não usar objetos dist.
mat2<- metaMDS(cem_lagos$Matrices, autotransform = T, try=200)
envfit(mat2, metadata[,c(4, 8:41, 52:55)])
plot(mat2, type="n")
points(mat2$points, col=metadata$catestadotrof)
points(mat2$points, col=metadata$CategBacia)
points(mat2$points, col=metadata$Cat_cla)
points(mat2$points, col=metadata$cat_pH)
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

# Rich
unlist(snoop[,6][snoop[,3]<=0.05])
# Shannon
unlist(snoop[,6][snoop[,4]<=0.05])
# Simp
unlist(snoop[,6][snoop[,5]<=0.05])
# Inv
unlist(snoop[,6][snoop[,6]<=0.05])



