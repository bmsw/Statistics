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
mite.h <- decostand (mite, "hellinger")
mite.xy.c <- scale(mite.xy, center = TRUE, scale = FALSE)

# Teste


distancias<-as.matrix(dist(mite.xy.c, "euclidian"))
d0<-distancias[upper.tri(distancias)]

fd<-function(x){
  h = 2*IQR(x)*(length(x)^(-(1/3)))
  bins=(max(x)-min(x))/h
  return(h)
}
# Distancia ideal por Friedman
fd(d0)
# Distancia ideal por Sturge
max(d0)/(1+log2(length(d0)))

#Scott
max(d0)/3.5*sd(d0)*length(d0)^-(1/3)

# Univariate spatial correlogram (based on Moran's I) =============

# Search for neighbours of all points within a radius of 0.7 m
# and multiples (i.e., 0 to 0.7 m, 0.7 to 1.4 m and so on).

plot.links(mite.xy, thresh = 0.46)
nb1 <- dnearneigh(as.matrix(mite.xy), 0, .46)
summary(nb1)

# Correlogram of substrate density
subs.dens <- mite.env[ ,1]
subs.correlog <- 
  sp.correlogram(nb1, 
                 subs.dens, 
                 order = 14, 
                 method = "I", 
	     zero.policy = TRUE)
print(subs.correlog, p.adj.method = "holm")

plot(subs.correlog)

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
  ylab("Moran's I")+xlab("Lags")+ggtitle("Correlogram", subtitle = NULL)
