library(sp)
library(rgdal)
library(raster)


# Create a grid of points within the bbox of the SpatialPolygonsDataFrame 
# colorado with decimal degrees as map units
grid <- makegrid(mapa, cellsize = 0.1)

# grid is a data.frame. To change it to a spatial data set we have to
grid <- SpatialPoints(grid, proj4string = CRS(proj4string(mapa)))
grid <- grid[mapa, ]

plot(mapa)
plot(grid, pch = ".", add = T)
