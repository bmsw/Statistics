source("/home/brunomsw/Documentos/DOUTORADO/flowDiv/flowDiv/flowDiv/flowDiv_internals.R")
#source("/home/brunomsw/Documentos/DOUTORADO/flowDiv/flowDiv/flowDiv/flowDiv_main.R")
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
setwd("/home/brunomsw/Documentos/DOUTORADO/DATA_FCS/Tara/Results")
# FASE 01: Importação e preparo dos datasets
#  Importar dataset OM_CompanionMODIFIED_ALL
OM_CompanionMODIFIED_ALL <- read.delim("~/Documentos/DOUTORADO/DATA_FCS/Tara/TARA METADATA/OM.CompanionMODIFIED_ALL.csv", sep="\t", dec=",", strip.white = T)




om.regexpr=regexpr("TARA_[0-9]{3}_[A-Z]{3}",OM_CompanionMODIFIED_ALL$`Sample.label..TARA_station._environmental.feature_size.fraction.`)
#ID's
om=regmatches(OM_CompanionMODIFIED_ALL$`Sample.label..TARA_station._environmental.feature_size.fraction.`,om.regexpr)
# Fractions
OM_fraction = regmatches(as.character(OM_CompanionMODIFIED_ALL$`Sample.label..TARA_station._environmental.feature_size.fraction.`),om.regexpr, invert = T)
fr=factor(gsub("^_", "", do.call( rbind, OM_fraction)[,2]))
# New DataSet
OM_CompanionMODIFIED_ALL=cbind(Station=om,Fractions=fr, OM_CompanionMODIFIED_ALL)


wksp=opc("../Bruno_Tara_ALL_Bact10-4.wsp")
######### flowDiv
tara.fd=flowDiv3(wksp, "Bact", do.plot = F, use.beads = T,ialpha = "shannon", beads = "Beads", static = F, transform = F, autotrans=F, psize = 1)

# By three
#tara_235=new("flowdiv", alfa=tara.fd$Alpha, pielou=tara.fd$Pielou,matrices=tara.fd$Matrices )
#tara_146=new("flowdiv", alfa=tara.fd$Alpha, pielou=tara.fd$Pielou,matrices=tara.fd$Matrices )


# By two
#tara_12=new("flowdiv", alfa=tara.fd$Alpha, pielou=tara.fd$Pielou,matrices=tara.fd$Matrices )
#tara_13=new("flowdiv", alfa=tara.fd$Alpha, pielou=tara.fd$Pielou,matrices=tara.fd$Matrices )
#tara_14=new("flowdiv", alfa=tara.fd$Alpha, pielou=tara.fd$Pielou,matrices=tara.fd$Matrices )
#tara_15=new("flowdiv", alfa=tara.fd$Alpha, pielou=tara.fd$Pielou,matrices=tara.fd$Matrices )
#tara_16=new("flowdiv", alfa=tara.fd$Alpha, pielou=tara.fd$Pielou,matrices=tara.fd$Matrices )
#tara_23=new("flowdiv", alfa=tara.fd$Alpha, pielou=tara.fd$Pielou,matrices=tara.fd$Matrices )
#tara_24=new("flowdiv", alfa=tara.fd$Alpha, pielou=tara.fd$Pielou,matrices=tara.fd$Matrices )
#tara_25=new("flowdiv", alfa=tara.fd$Alpha, pielou=tara.fd$Pielou,matrices=tara.fd$Matrices )
#tara_26=new("flowdiv", alfa=tara.fd$Alpha, pielou=tara.fd$Pielou,matrices=tara.fd$Matrices )
#tara_34=new("flowdiv", alfa=tara.fd$Alpha, pielou=tara.fd$Pielou,matrices=tara.fd$Matrices )
#tara_35=new("flowdiv", alfa=tara.fd$Alpha, pielou=tara.fd$Pielou,matrices=tara.fd$Matrices )
#tara_36=new("flowdiv", alfa=tara.fd$Alpha, pielou=tara.fd$Pielou,matrices=tara.fd$Matrices )
#tara_45=new("flowdiv", alfa=tara.fd$Alpha, pielou=tara.fd$Pielou,matrices=tara.fd$Matrices )
#tara_46=new("flowdiv", alfa=tara.fd$Alpha, pielou=tara.fd$Pielou,matrices=tara.fd$Matrices )
#tara_56=new("flowdiv", alfa=tara.fd$Alpha, pielou=tara.fd$Pielou,matrices=tara.fd$Matrices )


# By one
tara_1=new("flowdiv", alfa=tara.fd$Alpha, pielou=tara.fd$Pielou,matrices=tara.fd$Matrices )
tara_2=new("flowdiv", alfa=tara.fd$Alpha, pielou=tara.fd$Pielou,matrices=tara.fd$Matrices )
tara_3=new("flowdiv", alfa=tara.fd$Alpha, pielou=tara.fd$Pielou,matrices=tara.fd$Matrices )
tara_4=new("flowdiv", alfa=tara.fd$Alpha, pielou=tara.fd$Pielou,matrices=tara.fd$Matrices )
tara_5=new("flowdiv", alfa=tara.fd$Alpha, pielou=tara.fd$Pielou,matrices=tara.fd$Matrices )
tara_6=new("flowdiv", alfa=tara.fd$Alpha, pielou=tara.fd$Pielou,matrices=tara.fd$Matrices )



boxplot(list(tara_235@alfa, tara_146@alfa))

setClass("flowdiv",
         representation(
           alfa="numeric",
           pielou="numeric",
           matrices="matrix"
         )
)



# Preparing dataset
id.flowdiv_raw=names(tara.fd$Alpha)
id.flowdiv_raw_regexpr=regexpr("[0-9]{1,3} [a-z]{1,4}", id.flowdiv_raw, ignore.case = T)
id.flowdiv_clean=regmatches(id.flowdiv_raw,id.flowdiv_raw_regexpr, invert = F)
# Separates numbers and locals from flowDiv
## Numbers
fd.numbers=regexpr("[0-9]{1,3}", id.flowdiv_clean, ignore.case = T)
fd.numbers1 = regmatches(id.flowdiv_clean,fd.numbers)
# Adding "zeros" 
# One char
ones=nchar(fd.numbers1)==1
fd.numbers1[ones]=paste("00",sep="", fd.numbers1[ones])
# Two chars
twos=ones=nchar(fd.numbers1)==2
fd.numbers1[twos]=paste("0",sep="", fd.numbers1[twos])

#Fitting Pangea (tara1)
tara1=paste("TARA_", fd.numbers1, sep="") # Use tara1 for Panagea dataset
tara1=data.frame("Station"=tara1)
## Locals
fd.letters=regexpr("[a-z]{1,4}", id.flowdiv_clean, ignore.case = T)
fd.letters1= regmatches(id.flowdiv_clean,fd.letters)
# Fiting OM dataset (tara2)
fd.letters1[grep("SU", fd.letters1, ignore.case = T)]="SRF"
fd.letters1[grep("ME", fd.letters1, ignore.case = T)]="MES"
fd.letters1[grep("DCM", fd.letters1, ignore.case = T)]="DCM"
tara2=paste(as.matrix(tara1),"_", fd.letters1, sep="") 
tara2=data.frame("Station"=tara2)


########################################
# FASE 02: Analise dos dados
# Preparo da matriz final (de acordo com o dataset)
#cbind(rownames(tara.fd$Matrices), as.character(tara2$Station)) # Conferindo relação de linhas
rownames(tara.fd$Matrices)<-as.character(tara2$Station)
# Juntando datasets
OM2=join(tara2, OM_CompanionMODIFIED_ALL, by="Station")
# So fracoes .22-3 !!!! Legendas estão diferentes!(confira colunas de size fraction)
OM3=OM2[!is.na(OM2$Fractions),]
OM4=OM3[c(OM3$Fractions=="0.22-3"),]
#rownames(tara.fd$Matrices) %in% as.character(OM4$Station) # Checando linhas
#cbind(rownames(m1), as.character(OM4$Station)) # Checagem visual da organização
#Labels segundo as regioes
ocean.regions=factor(regmatches(OM4$`Ocean.and.sea.regions..IHO.General.Sea.Areas.1953...MRGID.registered.at.www.marineregions.com.`, regexpr("\\((.*)\\)", OM4$`Ocean.and.sea.regions..IHO.General.Sea.Areas.1953...MRGID.registered.at.www.marineregions.com.`)))
sample.regions=factor(regmatches(OM4$Station, regexpr("[A-Z]{3}$", OM4$Station)))
# Matriz TARA
m1=tara.fd$Matrices[rownames(tara.fd$Matrices) %in% as.character(OM4$Station),]

#Definindo os indices
rich=apply(m1, 1, function(x) specnumber(x, MARGIN = 1))
pielou=apply(m1, 1, function(x)diversity(x)/log(vegan::specnumber(x)))
shannon=apply(m1, 1, function(x) diversity(x,index = "shannon"))
simpson=apply(m1, 1, function(x) diversity(x,index = "simpson"))
invsimpson=apply(m1, 1, function(x) diversity(x,index = "invsimpson"))

OM5=cbind(OM4[,1:5],ocean.regions, sample.regions, "FC - Richness"=rich, "FC - Pielou"=pielou, "FC - Shannon"=shannon, "FC - Simpson"=simpson, "FC - invSimpson"=invsimpson, OM4[,18:52])
OM5=OM5[,complete.cases(t(OM5))] # Todas as colunas. Necessário?

# Fase 03: Plots

g1=ggplot(aes(y=`FC - Richness`, x=`sample.regions`), data=OM5)+geom_boxplot(aes(fill=`sample.regions`, alpha=.5))+ theme(legend.position="none")+ scale_x_discrete(name="")+geom_point()+geom_text(aes(label=OM5$Station.identifier..TARA_station..),hjust=-0.2, vjust=0, size=3)
g1
g2=ggplot(aes(y=`FC - Pielou`, x=`sample.regions`), data=OM5)+geom_boxplot(aes(fill=`sample.regions`, alpha=.5))+ theme(legend.position="none")+ scale_x_discrete(name="")+geom_point()+geom_text(aes(label=OM5$Station.identifier..TARA_station..),hjust=-0.2, vjust=0, size=3)
g2
g3=ggplot(aes(y=`FC - Shannon`, x=`sample.regions`), data=OM5)+geom_boxplot(aes(fill=`sample.regions`, alpha=.5))+ theme(legend.position="none")+ scale_x_discrete(name="")+geom_point()+geom_text(aes(label=OM5$Station.identifier..TARA_station..),hjust=-0.2, vjust=0, size=3)
g3
g4=ggplot(aes(y=`FC - Simpson`, x=`sample.regions`), data=OM5)+geom_boxplot(aes(fill=`sample.regions`, alpha=.5))+ theme(legend.position="none")+ scale_x_discrete(name="")+geom_point()+geom_text(aes(label=OM5$Station.identifier..TARA_station..),hjust=-0.2, vjust=0, size=3)
g4
g5=ggplot(aes(y=`FC - invSimpson`, x=`sample.regions`), data=OM5)+geom_boxplot(aes(fill=`sample.regions`, alpha=.5))+ theme(legend.position="none")+ scale_x_discrete(name="")+geom_point()+geom_text(aes(label=OM5$Station.identifier..TARA_station..),hjust=-0.2, vjust=0, size=3)
g5
g6=ggplot(aes(y=`FC...bacteria..cells.mL.`, x=`sample.regions`), data=OM5)+geom_boxplot(aes(fill=`sample.regions`, alpha=.5))+ theme(legend.position="none")+ scale_x_discrete(name="")+geom_point()+geom_text(aes(label=OM5$Station.identifier..TARA_station..),hjust=-0.2, vjust=0, size=3)
g6

g7=plot_grid(g1, g2, g3 ,g4 ,g5 ,g6)
g7

# Renomeando as colunas
OM6=OM5
newnames=c("St", "Fr", "Slb", "PgID", "StID", "OR", "SR", "FRc", "FPi", "FSh", "FSm","FInv", "mLat", "mLon", "mDep","mTemp","mOxi", "mDepS","mDepF","mDepN","mDep+O", "mDep-O", "miRich", "miDiv", "Chao", "Ace", "miShn","OGShn","OGRich", "OGEvn", "Hetr", "Aut","Bact","Pico", "mGT")
colnames(OM6)<-newnames
#Pairwise
scatterplotMatrix(OM6[,c(8:12, 23:34)], groups =OM5$sample.regions, by.groups=T, smooth=F, legend.plot = F, col=c("red", "green", "blue"))


# Correlation
cor.mtest <- function(mat, conf.level = 0.95){
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      tmp <- cor.test(mat[,i], mat[,j], conf.level = conf.level)
      p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
      lowCI.mat[i,j] <- lowCI.mat[j,i] <- tmp$conf.int[1]
      uppCI.mat[i,j] <- uppCI.mat[j,i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}

corr.test(OM6[,c(8:12, 23:34)])
res1 <- cor.mtest(OM6[,c(8:12, 23:34)], 0.95)
corrplot(cor(OM6[,c(8:12, 23:34)]), p.mat = res1[[1]], sig.level=.05, type = "lower", pch.col="red", method = "number", cl.cex=.8, number.cex = .7)

# PC
## TARA_082_SRF [71] seems an outlier! Removed from PCA
wine.pca <- prcomp(OM6[-71,c(8:35)], scale. = TRUE)

ggbiplot(wine.pca, obs.scale = 1, var.scale = 2,
         groups = OM6$SR[-71], ellipse = TRUE, circle = FALSE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')+
  scale_x_continuous(limits = c(-10, 10))+
  scale_y_continuous(limits = c(-10, 10))
########################################################################

# NMDS
# Plot from vegan
meta=metaMDS(m1[-71,], try=100, autotransform = F)
meta=metaMDS(decostand(m1[-71,], "hell"), try=100, autotransform = T)

plot(meta, type="n")
points(meta$points, col=as.numeric(sample.regions)+1)
fit=envfit(meta, OM6[-71,c(8:35)])
plot(fit, p.max=.05)
# plot from ggplot2
ggplot(aes(MDS1, MDS2, shape=ocean.regions[-71], colour=sample.regions[-71]), data=data.frame(meta$points))+geom_point(size=3)





