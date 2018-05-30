# script to analyze the functional diversity of coffee bird data
library(vegan)
library(cluster)
library(lattice)
source("/home/bruno/Downloads/FunctionalDiversity/FD_Calculators.R")
birdTraits <- read.csv("/home/bruno/Downloads/FunctionalDiversity/Bird_traits_annon.csv", header=TRUE)

head(birdTraits)
hist(birdTraits$Length_cm)
birdTraits$logLen=log2(birdTraits$Length_cm)
hist(birdTraits$logLen)


# lets drop a few of these because they are redundant or contain little info

birdTraitsr <- subset(birdTraits, select = c("Common", "logLen",
                                             "foodtype1", "forageHt",
                                             "sColor", "abun"))
# Now get abundance data from last week
birds <- read.csv("/home/bruno/Downloads/FunctionalDiversity/Coffee Bird List I.csv")
head(birds)
birds <- birds[, c(-1,-2)]
for (i in (2:length(birds))) {birds[,i]=ifelse(birds[,i]==0, NA,
                                               birds[,i])}
# must convert 0's to NA
birdsr=merge(birdTraitsr, birds, by="Common")
head(birdsr)
rownames(birdsr) <- birdsr[,"Common"];
head(birdsr)
# 
# # note on indexing. The brackets are a convenient way to index an
# object, array, or dataframe.
# # Today we will lump all our data in a single dataframe,
# “birdsr”, and use indexing to utilize parts of that dataframe.
# # The brackets contain two entries separated by a comma.
# # The first entry indexes the rows, and the second indexes the
# columns.
# # In you can use most any sort of expression in either spot. In
# the command above, we select only the column called “common” and
# assign those values to the row names. Below, we select just the
# farm abundance columns, or just the trait columns. When we work
# with convex hulls below, we will doubly index the dataframe
# birdsr: “birdsr[, trtcols][bighull, 2]”. See below for more.


FDAll=Calculate.FD(landfile=birdsr[,7:13],
                   traitfile=birdsr[,2:6])
QAll=Calculate.Q(landfile=birdsr[,7:13], traitfile=birdsr[,2:6])
FDAll
QAll
FDdf=data.frame(Site=names(FDAll), FD=as.numeric(FDAll))
Qdf=data.frame(Site=names(QAll), Q=as.numeric(QAll))
management <- read.csv("/home/bruno/Downloads/FunctionalDiversity/Coffee_Bird_Management.csv")
head(management)
management = merge(management, FDdf)
management = merge(management, Qdf)
management$ab=c('a','b','a','b','a','b','b')
quartz()
xyplot(FD~Management.index, groups=ab, data=management,
       auto.key=list(x=0.1, y=.95, border=TRUE))
summary(lm(FD~Management.index, data=management))
quartz()
xyplot(Q~Management.index, groups=ab, data=management,
       auto.key=list(x=0.1, y=.95, border=TRUE))
summary(lm(Q~Management.index, data=management))


# so lets plot the functional dendrograms for individual
#locations to see how FD compares between farms
# here is how we build the tree for all species present

cluster.distances <- daisy(
  birdsr[,2:6],
  metric="gower")
tree = hclust(cluster.distances)
tree <- as.dendrogram(tree)
quartz(width=6, height=12) # on windows use windows() instead of quartz()

par(mar = c(5, 4, 4, 18)) # this sets the margins of the plot as a number of lines, clockwise beginning at the bottom

plot(tree,
     main = "",
     horiz = TRUE,
     leaflab = "none", # this suppresses the automatic labeling of the leaves
)

# now we simply add leaf labels for the community (farm) of interest

mtext(as.character(row.names(birdsr)[order.dendrogram(tree)]), 4,
      las = 1, at = 1:length(birds$Common), # this orders the row names (=common names) as in the tree
      col = ifelse(is.na(birdsr$Forest_A)[order.dendrogram(tree)],
                   "grey50", "black")) # and this makes the leaf label grey if the species is missing, black if present let’s add a title too
title("functional dendrogram for Forest_A")
# 
# # you can now repeat this for each community (farm) you can open
# a new graphics window for each (quartz()) or write them to pdf
# pdf() try ?pdf for more info.
# # Question: how do the dendrograms of FD compare? Are there any
# clusters that are under or over represented on any farms?
#   # A third approach to quantifying functional diversity is to
#   calculate the convex hull volume that a community occupies in
# multidmensional trait space.
# # In this case we will calculate hull volume using just the 3
# numeric traits
# # here we will need the rgl and geometry packages
# # geometry will allow us to calculate the volume of n-dimensional
# convex hulls
# # rgl is a 3d viewer so we can see the results

library(rgl)
library(geometry)
library(reshape)
# first we need to change the format of birds (only the abundance data) from 'wide' to 'long'

birdsLong=melt(birdsr[,c(1,7:13)], id="Common",
               variable_name="community", na.rm=TRUE)
head(birdsLong)
rgl.open() # this opens a new rgl window

rgl.viewpoint(fov=10, zoom=1.5, theta=340, phi=0) # this sets the angle of view

rgl.bg(color="white") # this sets the background color
# 
# # so, here is the plan. First we will make a hull of all of the
# species, and compute its volume.
# # convex hull volume only works with continuous variables, so we
# will only use log(length), forage height, and abundance
# # Next we will make hulls of the individual communities and
# compare those volumes as a fraction of the entire community

head(birdsr)
# first we can contsruct and plot the full hull
trtcols=c("logLen","forageHt","abun")
rgl.clear( type = "shapes" )
# # First we use convhulln from the geometry package to find the
# convex hull of the coud of species points in 3d trait space.
# # The convex hull is simply the set of points that encloses all
# other points
# # That same hull can be defined as the set of facets or triangles
# (half-spaces in geometry-ese) that encloses the set of points.
# # In this case “bighull” is the set of triangles that represents
# the convex hull. Each column contains 3 numbers that represent
# the row number indices of the three points that bound each
# triangle in the original dataframe birdsr.

bighull = t(convhulln(birdsr[, trtcols], options="QJ"))
rgl.triangles(birdsr[, trtcols][bighull, 1], birdsr[,
                                                    trtcols][bighull, 2], birdsr[, trtcols][bighull,
                                                                                            3],col="blue",alpha=.8)
# # Here we are using rgl to plot the individual triangles that
# form the convex hull.
# # The three indices of birdsr (eg, birdsr[, trtcols][bighull, 1])
# each references one trait, such that the first 3 of each plot one
# facet, the second three the next facet, etc.


axis3d(c('x'), color= "black", labels=c(2,4,8,16,32))
axis3d(c('y'), color= "black", labels=c(1:5))
axis3d(c('z'), color= "black", labels=c(1:5))
title3d(main="Full community", sub=NULL, "Length (cm)", "Foraging
Height", "Abundance", color= "black", line=3)


# # if you click and drag of the RGL device (the hull) you can
# rotate it. Remind us to show you a movie made using this
# approach!
#   # now, we can plot each individual hull for the farms
#   # we will use a loop to make the code more compact
#   # this should make 7 windows with two 3d hulls in each.
#   # The light blue is the full community and the dark green is the
#   individual farm

vol=list()
for (i in levels(birdsLong$community)) { # i="Forest_A"
  rgl.open() # this opens a new rgl window
  rgl.viewpoint(fov=10, zoom=1.5, theta=340, phi=0) # this sets the angle of view
  rgl.bg(color="white") # this sets the background color
  rgl.triangles(birdsr[, trtcols][bighull, 1], birdsr[,
                                                      trtcols][bighull, 2], birdsr[, trtcols][bighull,
                                                                                              3],col="blue",alpha=.2)
  splist= birdsLong[birdsLong$community==i, "Common"]
  commhull = t(convhulln(birdsr[birdsr$Common%in%splist,
                                trtcols], options=paste("FS TO ",i,".txt", sep="")))
  rgl.triangles(birdsr[birdsr$Common%in%splist,
                       trtcols][commhull, 1], birdsr[birdsr$Common%in%splist,
                                                     trtcols][commhull, 2], birdsr[birdsr$Common%in%splist,
                                                                                   trtcols][commhull, 3],col="green",alpha=.8)
  axis3d(c('x'), color= "black", labels=c(2,4,8,16,32))
  axis3d(c('y'), color= "black", labels=c(1:5))
  axis3d(c('z'), color= "black", labels=c(1:5))
  title3d(main=i, sub=NULL, "Length (cm)", "Foraging Height",
          "Abundance", color= "black", line=3)
}

# well, that was fun. Now, notice that there are several rgl
# windows. You can rotate the hulls by clicking and dragging on
# them.
# But wait, what are the volumes of these hulls?
# Note that the volume is returned in the console, but there seems to be no way to get it as an object in R.
# Fortunately, there is a method to write the output to a file, which we can then read in as an R object.
management$vol=NULL # first make an empty container for the volumes

# now loop through the farms and calculate the volume for each.
for (i in levels(birdsLong$community))
{management[management$Site==i, "vol"]=scan(paste(i,".txt",
                                                  sep=""))[4]}
xyplot(vol~Management.index, groups=ab, data=management)
# Clearly, there is no effect of management intensity on this measure of diversity either.
# 
# # Further questions:
# # In your estimation, does farm management affect bird species or
# functional diversity?
#   # How might this analysis be improved?
#   # Are there any additional traits you would include if you could
#   obtain them?
  
