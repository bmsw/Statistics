library(ggplot2)
library(rgdal)
library(plyr)
library(rgeos)
library(ggmap)
library(sp)
library(raster)
library(fields)

# Mapa do RN (download do IBGE: https://mapas.ibge.gov.br/bases-e-referenciais/bases-cartograficas/malhas-digitais.html)
setwd("/home/bruno/Downloads/RN/")
# Lendo os arquivos 
mapa=readOGR(dsn=".","24MEE250GC_SIR", verbose = F)
plot(mapa)

# Create a grid of points within the bbox of the SpatialPolygonsDataFrame 
# colorado with decimal degrees as map units
grid <- makegrid(mapa, cellsize = 0.05)
# grid is a data.frame. To change it to a spatial data set we have to
grid <- SpatialPoints(grid, proj4string = CRS(proj4string(mapa)))
plot(grid)
#Subset
grid2 <- grid[mapa, ]
# Plot
plot(mapa)
plot(grid2, pch = ".", add = T)
#Krigagem

# Krigaem (exemplo com o pH)
kr2<- Krig(cbind(Lon.DECIMAL...., Lat.DECIMAL....), ph)
# Conferindo predições (não rodar)
# predict.Krig(kr2,cbind(Lon.DECIMAL...., Lat.DECIMAL....))
# Prevendo valores para o grid
pred<-predict.Krig(kr2,grid2@coords)
# Novo data.frame
new=data.frame(pred,lon=grid2@coords[,1], lat=grid2@coords[,2])


# Testes...
library(gridExtra)

plot1<-ggplot(aes(Lon.DECIMAL...., Lat.DECIMAL....), data=metadata) + geom_point(size=1) + coord_equal() + 
ggtitle("Points with measurements")

plot2<-ggplot(aes(lon, lat), data=new) + geom_point(size=1) + coord_equal() + 
ggtitle("Points at which to estimate")
# Similar to par=mfrow()
grid.arrange(plot1, plot2, ncol = 2)

ggplot(aes(x=lon, y=lat), data=new) + geom_tile(aes(fill=pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="red") +
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw()


