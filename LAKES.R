# ## Required Packages
library("betapart")
# library("reshape2")
library("vegan")
library("flowWorkspace")
# library("ggplot2")
# library("ggcyto")
# library("plyr")
library("car")
# library("corrplot")
# library("pangaear")
# #library("ggbiplot")
# library("psych")
# library("cowplot")
# library("RColorBrewer")
# library("dunn.test")
# library("seriation")
# library("gvlma")
# library("RVAideMemoire")
# library("gplots")
# library("Imap") # Para a função GeoDistance
library(scatterpie)
# Sourcing
source("/home/bruno/Documentos/DOC/flowDiv_SUBMISSION/Packages_and_Codes/flowDiv2/R/flowDiv_internalsCORRECTED.R")
source("/home/bruno/Documentos/DOC/flowDiv_SUBMISSION/Packages_and_Codes/flowDiv2/R/flowDiv_main.R")


## Diretorio de trabalho e Metadados
setwd("/home/bruno/Documentos/DOC/MAIN/100LAGOS/FCS/")

## Nao usar blobo abaixo - novas variaveis salvas NEW_ID e NEW_MF

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

## Nao usar blobo acima - novas variaveis salvas NEW_ID e NEW_MF


metadata=read.csv("NEW_MF.csv")
ids=read.csv("NEW_ID_DILUTION.csv")

# flowDiv com as diluições
cem_lagos = flowDiv("Bruno_100Lagos_Beadsv10.wsp", "Bact", dilutions = ids$DILUTION)
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

# Krigagem

library(rgdal)
library(plyr)
library(rgeos)
library(ggmap)
library(sp)
library(raster)
library(fields)
library(gridExtra)

# Mapa do RN (download do IBGE: https://mapas.ibge.gov.br/bases-e-referenciais/bases-cartograficas/malhas-digitais.html)
setwd("/home/bruno/Downloads/RN/")
# Lendo os arquivos 
mapa=readOGR(dsn=".","24MEE250GC_SIR", verbose = F)
plot(mapa)
# Create a grid of points within the bbox of the SpatialPolygonsDataFrame 
# colorado with decimal degrees as map units
grid <- makegrid(mapa, cellsize = 0.01)
# grid is a data.frame. To change it to a spatial data set we have to
grid <- SpatialPoints(grid, proj4string = CRS(proj4string(mapa)))
plot(grid)
#Subset
grid2 <- grid[mapa, ]
# Plot
plot(mapa)
plot(grid2, pch = ".", add = T)
#Krigagem

kr_rich<- with(dataall, Krig(cbind(Lon, Lat), rich))
kr_pielou<- with(dataall, Krig(cbind(Lon, Lat), pielou))
kr_shannon<- with(dataall, Krig(cbind(Lon, Lat), shannon))

# Prevendo valores para o grid
pred_rich<-predict.Krig(kr_rich,grid2@coords)
pred_pielou<-predict.Krig(kr_pielou,grid2@coords)
pred_shannon<-predict.Krig(kr_shannon,grid2@coords)

# Novo data.frame
map_kr=data.frame(pred_rich,pred_pielou,pred_shannon,lon=grid2@coords[,1], lat=grid2@coords[,2])


# Rich
plotr<-ggplot(aes(x=lon, y=lat), data=map_kr) + geom_tile(aes(fill=pred_rich)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="red", name="Richness") +
  scale_x_continuous() + scale_y_continuous() +
  labs(x="Longitude", y="Latitude")+theme_bw()

# Pielou
plotp<-ggplot(aes(x=lon, y=lat), data=map_kr) + geom_tile(aes(fill=pred_pielou)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="red", name="Pielou Index") +
  scale_x_continuous() + scale_y_continuous() +
  labs(x="Longitude", y="Latitude")+theme_bw()

# Shannon
plots<-ggplot(aes(x=lon, y=lat), data=map_kr) + geom_tile(aes(fill=pred_shannon)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="red", name="Shannon Index") +
  scale_x_continuous() + scale_y_continuous() +
  labs(x="Longitude", y="Latitude")+theme_bw()

#Plot

grid.arrange(plotr, plotp, plots, ncol = 3)


######################
library("ggplot2")  
library("ggrepel")  
library("maptools")  
library("broom")

mapa.tidy <- tidy(mapa) 
mapa.tidy <- merge(mapa.tidy, mapa@data)

g1<-ggplot(mapa.tidy) +  
  aes(long, lat, group = group) +
  geom_polygon(fill="gray70") +
  coord_equal()+
  labs(x="Longitude", y="Latitude")+theme_bw()

g1+geom_scatterpie(aes(x=Lon, y=Lat),data=dataall, cols=c("N.P", "a430"))

######################

# Particionamento beta
bpart=bray.part(m1)


nest<-as.matrix(bpart[[1]])
turn<-as.matrix(bpart[[2]])
bray<-as.matrix(bpart[[3]])

data_bray<-data.frame(Nestdness=nest[1,], Turnover=turn[1,], Lon=dataall$Lon, Bray=bray[1,], Lat=dataall$Lat)
data_bray2<-data_bray[-1,]

#g1+geom_scatterpie(aes(x=Lon, y=Lat), data=data_bray2, cols=c("Nestdness", "Turnover"))+theme(legend.title = element_blank())
# Raio com escala
# g1+geom_scatterpie(aes(x=Lon, y=Lat, r=Bray/14), data=data_bray2, cols=c("Nestdness", "Turnover"))+theme(legend.title = element_blank())+
#   geom_scatterpie_legend(data_bray2$Bray/14 , x=-39, y=-7, n=2, labeller=function(x) round(14*x, digits = 1))
# 



#g1+geom_point(aes(x=Lon, y=Lat), data=data_bray2, inherit.aes = FALSE)

#

# data_bray<-data.frame(Nestdness=nest[1,], Turnover=turn[1,], Lon=dataall$Lon, Lat=dataall$Lat)
# data_bray2<-data_bray[-1,]
# 
# g1+geom_scatterpie(aes(x=Lon, y=Lat), data=data_bray2, cols=c("Nestdness", "Turnover"))+
#   geom_point(aes(x=Lon, y=Lat), data=data_bray[1,], inherit.aes = FALSE, size=9, pch="\u2605", colour="yellow1")+
#   theme(legend.title = element_blank())
# 


for (i in 1:5){
  data_bray<-data.frame(Nestdness=nest[i,], Turnover=turn[i,], Lon=dataall$Lon, Bray=bray[i,], Lat=dataall$Lat)
  
  data_bray2<-data_bray[-i,]
  plot1<-g1+geom_scatterpie(aes(x=Lon, y=Lat, r=Bray/14), data=data_bray2, cols=c("Nestdness", "Turnover"))+
    geom_point(aes(x=Lon, y=Lat), data=data_bray[i,], inherit.aes = FALSE, size=9, pch="\u2605", colour="yellow1")+
    theme(legend.title = element_blank())+
    ggtitle(paste("Lago: ", metadata$name[i]))+
    geom_scatterpie_legend(data_bray2$Bray/14 , x=-39, y=-7, n=2, labeller=function(x) round(14*x, digits = 1))
  
  plot(plot1)
  
}







# Nova variavel
#names=c("a","d","d","e","f",colnames(metadata))
######### ALFA
mydata=dataall # New dataset (to not override!)
# Using justo 10:20 >>> Change it!

for (i in 10:20){
  bc=boxCox(lm(mydata[,i]~mydata$pielou), family="yjPower")
  best=with(bc, x[which.max(y)])
  mydata[,ncol(mydata)+1]<-yjPower(mydata[,i], best)
  colnames(mydata)[ncol(mydata)]<-paste("TransP", colnames(mydata)[i], sep="_")
}

##################

## Função para verificar a adequação dos modelos lineares (gvlma OK pra todos os testes)
gvbruno<-function(a){
  soma=sum(a$GlobalTest$GlobalStat4$Decision, a$GlobalTest$DirectionalStat1$Decision, a$GlobalTest$DirectionalStat2$Decision, a$GlobalTest$DirectionalStat3$Decision, a$GlobalTest$DirectionalStat4$Decision)
  return(soma==0)}

# Verificando gvlma
gr=c(); gs=c(); gp=c(); gsi=c(); ginv=c()
for(i in 8:47) gr[i-7]=gvbruno(gvlma(lm(rich~metadata[,i])))
for(i in 8:47) gp[i-7]=gvbruno(gvlma(lm(pielou~metadata[,i])))
for(i in 8:47) gs[i-7]=gvbruno(gvlma(lm(shannon~metadata[,i])))
for(i in 8:47) gsi[i-7]=gvbruno(gvlma(lm(simpson~metadata[,i])))
for(i in 8:47) ginv[i-7]=gvbruno(gvlma(lm(invsimpson~metadata[,i])))


# Verificando Correlações

r=c(); p=c(); sh=c(); si=c(); inv=c()
for(i in 8:47) r[i-7]=cor.test(rich, metadata[,i])$p.value
for(i in 8:47) p[i-7]=cor.test(pielou, metadata[,i])$p.value
for(i in 8:47) sh[i-7]=cor.test(shannon, metadata[,i])$p.value
for(i in 8:47) si[i-7]=cor.test(simpson, metadata[,i])$p.value
for(i in 8:47) inv[i-7]=cor.test(invsimpson, metadata[,i])$p.value

# Dados das correlações <= 0.05 e gvlma's correspondentes
correls=cbind(cbind(p, r, sh, si, inv)<=0.05, gp, gr, gs,gsi, ginv, names[8:47])
correls<-data.frame(correls)
# Organizando tipos de dados para booleanos
##
correls$p<-as.logical(correls$p)
correls$r<-as.logical(correls$r)
correls$sh<-as.logical(correls$sh)
correls$inv<-as.logical(correls$inv)
correls$si<-as.logical(correls$si)
##
correls$gp<-as.logical(correls$gp)
correls$gr<-as.logical(correls$gr)
correls$gs<-as.logical(correls$gs)
correls$gsi<-as.logical(correls$gsi)
correls$ginv<-as.logical(correls$ginv)

# Snooping: Verificando modelos significativos && válidos (segundo gvlma)
snoop1=cbind(p=(correls$p)&(correls$gp), r=(correls$r)&(correls$gr), sh=(correls$sh)&(correls$gs), si=(correls$si)&(correls$gsi), inv=(correls$inv)&(correls$ginv))
snoop2=cbind(snoop1, names[8:47])
##
snoop3=snoop2[!apply(snoop1,1, sum)==0,]
# Análises visuais:
par(mfrow=c(3,3))
for (i in 1:14)plot(dataall$pielou, dataall[,snoop3[i,6]], ylab=snoop3[i,6])

# Texto
par(mfrow=c(3,3))
for (i in 1:14){plot(dataall$pielou, dataall[,snoop3[i,6]], ylab=snoop3[i,6], type="n"); text(dataall$pielou, dataall[,snoop3[i,6]], labels =metadata$lakeID)}

# regressões (Pielou): pH, Temp e LogCLA


ggplot(aes(pielou, ph), data=dataall)+geom_point()+geom_smooth(method='lm')+
  xlab("Pielou Index")+ylab("ph")
par(mfrow=c(2,2))
plot(lm(ph~pielou, dataall), main="pH")


ggplot(aes(pielou, temp.ºC), data=dataall)+geom_point()+geom_smooth(method='lm')+
  xlab("Pielou Index")+ylab("Temperature") # Not good
par(mfrow=c(2,2))
plot(lm(temp.ºC~pielou, dataall), main="Temperature")


ggplot(aes(pielou, LOGCLA.1), data=dataall)+geom_point()+geom_smooth(method='lm')+
  xlab("Pielou Index")+ylab("Log10 Chla")

par(mfrow=c(2,2))
plot(lm(LOGCLA.1~pielou, dataall), main="Chla")


######### BETA
bray=vegdist(m1)
m2=data.frame(name=metadata$name, lat=metadata$Lat.DECIMAL...., lon=metadata$Lon.DECIMAL....)
##########################################
ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
  # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
  # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.
  
  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}

GeoDistanceInMetresMatrix <- function(df.geopoints){
  # Returns a matrix (M) of distances between geographic points.
  # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
  # (df.geopoints$lat[j], df.geopoints$lon[j]).
  # The row and column names are given by df.geopoints$name.
  
  GeoDistanceInMetres <- function(g1, g2){
    # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
    # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
    # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
    # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
    # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
    DistM <- function(g1, g2){
      require("Imap")
      return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="m")))
    }
    return(mapply(DistM, g1, g2))
  }
  
  n.geopoints <- nrow(df.geopoints)
  
  # The index column is used to ensure we only do calculations for the upper triangle of points
  df.geopoints$index <- 1:n.geopoints
  
  # Create a list of lists
  list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})
  
  # Get a matrix of distances (in metres)
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")
  
  # Set the row and column names
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name
  
  return(mat.distances)
}
##########################################
# Matriz de distanca geografica
geo=round(GeoDistanceInMetresMatrix(m2) / 1000)
# Particionamento beta
bpart=bray.part(m1)
# Mantel
mantel(bpart[[1]], as.dist(geo))
mantel(bpart[[2]], as.dist(geo)) # =)
mantel(bpart[[3]], as.dist(geo))
