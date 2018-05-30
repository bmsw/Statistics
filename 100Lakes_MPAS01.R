# Maps

load("~/Data100Lakes.RData")


# Krigagem

library(rgdal)
library(plyr)
library(rgeos)
library(ggmap)
library(sp)
library(raster)
library(fields)
library(gridExtra)
library(ape)
library(spdep)
library(ade4)
library(adegraphics)
library(adespatial)
library(vegan)
library(SoDA) # for geoXY()

# Mapa do RN (download do IBGE: https://mapas.ibge.gov.br/bases-e-referenciais/bases-cartograficas/malhas-digitais.html)
setwd("/home/bruno/Downloads/RN/")
# Lendo os arquivos 
mapa=readOGR(dsn=".","24MEE250GC_SIR", verbose = F)
#plot(mapa)
# Create a grid of points within the bbox of the SpatialPolygonsDataFrame 
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
# Scatterpies
library("ggplot2")  
library("ggrepel")  
library("maptools")  
library("broom")
library("scatterpie")
library("betapart")

mapa.tidy <- tidy(mapa) 
mapa.tidy <- merge(mapa.tidy, mapa@data)

# Particionamento beta
bpart=bray.part(m1)

# Nes dataframe
nest<-as.matrix(bpart[[1]])
turn<-as.matrix(bpart[[2]])
bray<-as.matrix(bpart[[3]])

# Change n!
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

# Mantel correlogram
# Get cartesian coordinates from Lat and Lon

dataall.cart<-geoXY(dataall$Lat, dataall$Lon, unit=1000)
# Distanc class = 25km
mantcor<-function(ind, titulo){
  (lakes.correlog <- 
     mantel.correlog(bpart[[ind]], r.type = "spearman",
                     XY = dataall.cart, 
                     nperm = 999, break.pts = seq(0, 300, 25), cutoff = T))
  print(summary(lakes.correlog))
  
  # Number of classes
  print(lakes.correlog$n.class) # or: mite.correlog[2]
  # Break points
  print(lakes.correlog$break.pts) # or: mite.correlog[3]
  
  plot(lakes.correlog)
  title(titulo)
  
  # newd<-data.frame(lakes.correlog$mantel.res[1:4,])
  # sigs<-newd$Pr.corrected.[1:4]<=.05
  # shapes<-ifelse(sigs==TRUE, 8,16)
  # 
  # p1<-ggplot(aes(x=class.index, y=Mantel.cor), data=newd)+
  #   geom_line() +
  #   geom_point(shape=shapes, size=4) +
  #   geom_hline(yintercept=0, colour="red", size=.5)+
  #   scale_x_continuous(breaks=seq(0, 150, 50))+
  #   theme_bw()+
  #   ylab("Mantel Correlation")+xlab("Distance Class Index")+ggtitle("Mantel Correlogram", subtitle = titulo)+
  #   theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  # plot(p1)
  
}


mantcor(3, "Bray-Curtis")
mantcor(2, "Nestedness")
mantcor(1, "Turnover")

# Alphas

# Teste


distancias<-as.matrix(dist(dataall.cart, "euclidian"))
d0<-distancias[upper.tri(distancias)]
# Friedman-Diaconis
fd<-function(x){
  h = 2*IQR(x)*(length(x)^(-(1/3)))
  bins=(max(x)-min(x))/h
  return(bins)
}
# Distancia ideal por Friedman
fd(d0)
# Distancia ideal por Sturge
(max(d0)-min(d0))/(1+log2(length(d0)))
#Scott
(max(d0)-min(d0))/(3.5*sd(d0)*length(d0)^-(1/3))

# Sugestão: 25km

# Univariate spatial correlogram (based on Moran's I) =============

#25 Km
source("/home/bruno/Documentos/DOC/R_Scripts_Numerical_Ecology/NEwR-2ed_code_data/NEwR2-Functions/plot.links.R")
plot.links(dataall.cart, thresh = 25)
nb1 <- dnearneigh(as.matrix(dataall.cart), 0, 25)
summary(nb1)

# Correlogram 
correlg<-function(ind, ord=8){
  subs.dens <- dataall[,ind]
  subs.correlog <- sp.correlogram(nb1, subs.dens,order = ord , method = "I",zero.policy = TRUE)
  sigs<-print(subs.correlog, p.adj.method = "holm")[,5]
  plot(subs.correlog, main="Correlogram")
  #ggplot
  subs.correlog.df<-data.frame(lags=seq(1,nrow(subs.correlog$res)),estimate=subs.correlog$res[,1], signif=sigs)
    ggplot(aes(lags, estimate), data=subs.correlog.df)+
    geom_line()+
    geom_point(size=2, colour=ifelse(sigs<=0.05, "red", "black"))+
    scale_x_continuous(breaks = c(1:nrow(subs.correlog.df)))+
    scale_y_continuous(breaks = round(seq(-.8, .5, by = 0.05),1))+
    geom_hline(yintercept=-1/(nrow(distancias)-1), colour="black", size=.5, lty=1)+
    theme_bw()+
    ylab("Moran's I")+xlab("Lags")+ggtitle("Correlogram", subtitle = NULL)+
    theme(plot.title = element_text(hjust = 0.5))+
    annotate("text", x = 7.5, y = .2, label = "d = 25Km", size=3)
  
  
}
correlg("shannon")#
correlg("pielou")
correlg("rich")#

