# Diversidade citométrica 100 Lagos RN
## Script 04/?? : Functional
####Pacotes utilizados
library("vegan")
library("flowWorkspace")
library("betapart")
library("car")
library("psych")
####################### Códigos-fonte ##################################
source("/home/bruno/Documentos/DOC/flowDiv_SUBMISSION/Packages_and_Codes/flowDiv2/R/flowDiv_internalsCORRECTED.R")
source("/home/bruno/Documentos/DOC/flowDiv_SUBMISSION/Packages_and_Codes/flowDiv2/R/flowDiv_main2_2018teste.R")


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
2

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



library(FD)
traits<- cem_lagos$ids
abund<-m1
# Functional indices
teste=dbFD(traits, abund)
#enviromental
envf<-dataall[,c(9,13:37,48:50)]
colnames(envf)[c(8, 2, 28, 26, 9, 1, 29, 5)]=c("Abd", "RP", "Chla", "N.P", "pH", "Prec", "CO2", "PB")


nmds0=cca(abund~1,envf, scale=T)
nmds1=cca(abund~.,envf, scale=T)
best=ordistep(nmds0, nmds1)
plot(best, scaling=1)
anova(best, by = "terms")

library(vegan3d)
# Manual plot
plot(best, scaling=1, type="n")
# All
text(scores(best, choices = c(1,2), display="species", scaling=1))
#Chla
points(scores(best, choices = c(1,2), display="species", scaling=1), pch=15, col=cem_lagos$ids$PerCP.Cy5.5.H)
# Size
points(scores(best, choices = c(1,2), display="species", scaling=1), pch=15, col=cem_lagos$ids$SSC.H)
# DNA
points(scores(best, choices = c(1,2), display="species", scaling=1), pch=15, col=cem_lagos$ids$FITC.H)
#Sites
#points(scores(best, choices = c(1,2), display="sites", scaling=1))

bip <- scores(best, choices = 1:2, display = "bp", scaling=1)
escala=ordiArrowMul(best, display = "bp", scaling=1)
bip.scl=bip*escala
labs <- rownames(bip)
(bip.lab <- ordiArrowTextXY(bip.scl, rescale = FALSE, labels = labs))
arrows(0, 0, bip.scl[,1], bip.scl[,2], length = 0.1)
text(bip.lab, labels = labs)
title("Chla")

# RAW PIPE


myworkspaces="Bruno_100Lagos_Beadsv10.wsp"
gate.name="Bact"
do.plot = T
static = F
dilutions=ids$DILUTION
transform = T
autotrans=T
psize = 1.08
ialpha="invsimpson"; ibeta="bray"; pmax=0.05; use.beads=FALSE; beads=NULL; transform.hell=FALSE; dimension=20

classes=sapply(myworkspaces, function(x) return(class(x)), simplify = T)
gsets=which(classes=="GatingSet")
if(length(gsets)>0){
  wksp1=myworkspaces[-gsets]
  gtsets=myworkspaces[gsets]
  wksp2<-opc(wksp1)
  wksp<-append(wksp2, gtsets)

} else (wksp=opc(myworkspaces))
nnn=nn(wksp, nod=gate.name,use.beads=use.beads, nod2=beads)



ops<-unique(unlist(lapply(unlist(nnn$nodesample), colnames)))
#selection<-select.list(ops, c(1:length(ops)), multiple = T, title="Please select channels for use")
selection=c("SSC-H","FITC-H","PerCP-Cy5-5-H")
fixed=NULL
binss=lapply(nnn$nodesample, function(x) lapply(x, function(y)lapply(selection, function(x)Freedman.Diaconis(exprs(y)[,x]))))
suggested.bins=round(median(unlist(binss)))

message(paste("Suggested number of bins:", suggested.bins, "\n"), "How many bins do you want to use?")
#nbins=scan(nmax=1, quiet=T)
nbins=2
myseqs <- seq_fun(nnn$nodesample, selection, nbins, fixed)

myseqs

ratio=order(envf$PB)
library(ggplot2)
gs<-list()
for (i in 91:94){
  tubed=data.frame(exprs(nnn$nodesample[[1]][[ratio[i]]]))
  g1=ggplot(tubed)+
    aes_string(x="SSC.H", y="FITC.H")+
    geom_hex(bins = 100, na.rm = T) +
    scale_fill_gradientn("", colours = rev(gray.colors(3, end = 4/6)))+
    scale_y_continuous(limits= c(myseqs$max.min.table["FITC-H", "minimum"], myseqs$max.min.table["FITC-H", "maximum"])) +
    scale_x_continuous(limits= c(myseqs$max.min.table["SSC-H", "minimum"], myseqs$max.min.table["SSC-H", "maximum"])) +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line = element_blank(),
          axis.title=element_text(size=20),
          # axis.line = element_line(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    theme(legend.position="none")+
    geom_hline(yintercept=c(-0.144616357), linetype="dashed", color = "red")+
    geom_vline(xintercept=c(-0.25785841), linetype="dashed", color = "red")
  #geom_point(data=clusters2, aes_string(y="FL1.H", x="SSC.H"), colour=clusters2$clust, alpha=.5, shape=15, size=1.83)
  gs[[i-90]]<-g1
}

marrangeGrob(gs, nrow=2, ncol=2)
ggsave(marrangeGrob(gs, nrow=2, ncol=2),filename = "Higher", device = "pdf")
