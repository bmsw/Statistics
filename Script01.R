# Diversidade citométrica 100 Lagos RN
## Script 01/?? : Importação de metadados, geração de índices de diversidade
## e correlações
####Pacotes utilizados
library("vegan")
library("flowWorkspace")
library("betapart")
library("car")
library("psych")
####################### Códigos-fonte ##################################
source("/home/bruno/Documentos/DOC/flowDiv_SUBMISSION/Packages_and_Codes/flowDiv2/R/flowDiv_internalsCORRECTED.R")
source("/home/bruno/Documentos/DOC/flowDiv_SUBMISSION/Packages_and_Codes/flowDiv2/R/flowDiv_main.R")


## Diretorio de trabalho e Metadados
setwd("/home/bruno/Documentos/DOC/MAIN/100LAGOS/FCS/")

## Nao usar bloco abaixo - novas variaveis salvas NEW_ID e NEW_MF

# ids=read.csv("../METADATA/ID_DILUTION.csv", sep="\t")
# metadata=read.csv("../METADATA/MF.csv", dec=",")
# metadata=metadata[,-c(26:57)]

# Primeira: registrar ID's para organização das linhas dos metadados
# cem_lagos = flowDiv("Bruno_100Lagos_Beadsv10.wsp", gate.name = "Bact", beads = "Beads", use.beads = T )
#
# # Limpeza das planilhas: Alguns lagos não constam nas planilhas
# lakes<-unlist(lapply(names(cem_lagos$Alpha), function(x)substr(x, 1, 23)))
# ids=ids[ids$FCS%in%lakes,]
# metadata=metadata[metadata$lakeID%in%ids$LAGO,]
# # Reorganizando
# match(metadata$lakeID, ids$LAGO)
# metadata=metadata[match(ids$LAGO, metadata$lakeID),]
# metadata=metadata[,complete.cases(t(metadata))]
# write.csv(ids,"NEW_ID_DILUTION.csv", row.names = F)
# write.csv(metadata,"NEW_MF.csv", row.names = F)

## Nao usar bloco acima - novas variaveis salvas NEW_ID e NEW_MF


metadata=read.csv("NEW_MF.csv")
ids=read.csv("NEW_ID_DILUTION.csv")

# flowDiv com as diluições
cem_lagos = flowDiv("Bruno_100Lagos_Beadsv10.wsp", "Bact", dilutions = ids$DILUTION)
2 3 5
43

# Removendo Outlier: Riacho da Cruz (94)
metadata=metadata[-94,]
m1=cem_lagos$Matrices[-94,]
ids=ids[-94,]

# Clean 2
metadata<-metadata[,-c(6,7, 19:25)]
colnames(metadata)[6:7]<-c("Lat", "Lon")

# Definindo novos indices
rich=apply(m1, 1, function(x) specnumber(x, MARGIN = 1))
pielou=apply(m1, 1, function(x)diversity(x)/log(vegan::specnumber(x)))
shannon=apply(m1, 1, function(x) diversity(x,index = "shannon"))
simpson=apply(m1, 1, function(x) diversity(x,index = "simpson"))
invsimpson=apply(m1, 1, function(x) diversity(x,index = "invsimpson"))

# Novo dataframe
dataall<-cbind(rich, pielou, shannon, simpson, invsimpson, metadata)
dataall$Lat=-1*dataall$Lat
dataall$Lon=-1*dataall$Lon

# Correlações
cors<-corr.test(data.matrix(dataall), method = "spearman", use = "complete.obs", adjust = "none")
ps<-cors$p # p values
rs<-cors$r # statistics
rs[ps>0.05]<-NA # substituição por NA (apenas para faciliatr a analise visual)
rs

