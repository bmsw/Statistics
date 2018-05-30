### CHAPTER 7: SPATIAL ANALYSIS
###
### Online supporting material for: 
### Borcard D., Gillet F. & Legendre P. Numerical Ecology with R, 
### Second Edition, Springer, 2018

# Load packages, functions and data ===============================
library(ape)
library(spdep)
library(ade4)
library(adegraphics)
library(adespatial)
library(vegan)
library(SoDA) # for geoXY()


# Source additional functions that will be used later in this
# Chapter. Our scripts assume that files to be read are in
# the working directory.
source("/home/bruno/Documentos/DOC/R_Scripts_Numerical_Ecology/NEwR-2ed_code_data/NEwR2-Functions/plot.links.R")
source("/home/bruno/Documentos/DOC/R_Scripts_Numerical_Ecology/NEwR-2ed_code_data/NEwR2-Functions/sr.value.R")
source("/home/bruno/Documentos/DOC/R_Scripts_Numerical_Ecology/NEwR-2ed_code_data/NEwR2-Functions/quickMEM.R")
source("/home/bruno/Documentos/DOC/R_Scripts_Numerical_Ecology/NEwR-2ed_code_data/NEwR2-Functions/scalog.R")

# Load the oribatid mite data. The file mite.Rdata is assumed 
# to be in the working directory.
load("/home/bruno/Documentos/DOC/R_Scripts_Numerical_Ecology/NEwR-2ed_code_data/NEwR2-Data/mite.RData")

# Transform the data
m1.h <- decostand (m1, "hellinger")

m2=data.frame(name=dataall$lakeID, lat=dataall$Lat, lon=dataall$Lon)

geo=round(GeoDistanceInMetresMatrix(m2) / 1000)



# Get cartesian coordinates from Lat and Lon
dataall.cart<-geoXY(dataall$Lat, dataall$Lon, unit=1000)
dataall.xy.c <- scale(dataall.cart, center = TRUE, scale = FALSE)
dist2<-vegdist(dataall.xy.c, "euclidian")

# Teste


distancias<-as.matrix(dist(dataall.xy.c, "euclidian"))
d0<-distancias[upper.tri(distancias)]
# Friedman-Diaconis
fd<-function(x){
  h = 2*IQR(x)*(length(x)^(-(1/3)))
  bins=(max(x)-min(x))/h
  return(h)
}
# Distancia ideal por Friedman
fd(d0)
# Distancia ideal por Sturge
(max(d0)-min(d0))/(1+log2(length(d0)))

#Scott
max(d0)/(3.5*sd(d0)*length(d0)^-(1/3))
# Sugestão: 20km?

# Univariate spatial correlogram (based on Moran's I) =============

# Search for neighbours of all points within a radius of 0.7 m
# and multiples (i.e., 0 to 0.7 m, 0.7 to 1.4 m and so on).
#25 Km
plot.links(dataall.cart, thresh = 20.8)
nb1 <- dnearneigh(as.matrix(dataall.cart), 0, 20.8)
summary(nb1)

# Correlogram 
correlg<-function(ind, ord=8){
  subs.dens <- dataall[,ind]
  subs.correlog <- sp.correlogram(nb1, subs.dens,order = ord , method = "I",zero.policy = TRUE)
  print(subs.correlog, p.adj.method = "holm")
  plot(subs.correlog)
}
correlg("shannon", 10)#
correlg("pielou", 10)
correlg("rich", 10)#
# ggplot

subs.correlog.df<-data.frame(lags=seq(1,nrow(subs.correlog$res)),estimate=subs.correlog$res[,1], variance=2*sqrt(subs.correlog$res[,3]))


pd <- position_dodge(0.1) # move them .05 to the left and right
ggplot(aes(lags, estimate), data=subs.correlog.df)+
  geom_errorbar(aes(ymin=estimate-variance, ymax=estimate+variance), width=.3, position=pd)+
  geom_point(size=2)+
  scale_x_continuous(breaks = c(1:nrow(subs.correlog.df)))+
  scale_y_continuous(breaks = round(seq(-.8, .5, by = 0.05),1))+
  geom_hline(yintercept=-1/(nrow(distancias)-1), colour="red", size=1)+
  theme_bw()+
  ylab("Moran's I")+xlab("Lags")+ggtitle("Correlogram", subtitle = NULL)+
  theme(plot.title = element_text(hjust = 0.5))



# Mantel correlogram ====================

mite.h.D1 <- bray.part(m1)

mantcor<-function(ind, titulo){
(mite.correlog <- 
    mantel.correlog(mite.h.D1[[ind]], r.type = "spearman",
                    XY = dataall.cart, 
                    nperm = 999, break.pts = seq(0, 300, 50), cutoff = T))
summary(mite.correlog)

# Number of classes
mite.correlog$n.class # or: mite.correlog[2]
# Break points
mite.correlog$break.pts # or: mite.correlog[3]

plot(mite.correlog)
title(titulo)

newd<-data.frame(mite.correlog$mantel.res[1:4,])
sigs<-newd$Pr.corrected.[1:4]<=.05
shapes<-ifelse(sigs==TRUE, 8,16)

p1<-ggplot(aes(x=class.index, y=Mantel.cor), data=newd)+
  geom_line() +
  geom_point(shape=shapes, size=4) +
  geom_hline(yintercept=0, colour="red", size=.5)+
  scale_x_continuous(breaks=seq(0, 150, 50))+
  theme_bw()+
  ylab("Mantel Correlation")+xlab("Distance Class Index")+ggtitle("Mantel Correlogram", subtitle = titulo)+
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
plot(p1)

}


mantcor(3, "Bray-Curtis")
mantcor(2, "Nestedness")
mantcor(1, "Turnover")


