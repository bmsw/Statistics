# Load packages, functions and data ===============================
library(ape)
library(spdep)
library(ade4)
library(adegraphics)
library(adespatial)
library(vegan)
library(SoDA) # for geoXY()


# Source additional functions that will be used later 
source("/home/bruno/Documentos/DOC/R_Scripts_Numerical_Ecology/NEwR-2ed_code_data/NEwR2-Functions/plot.links.R")
source("/home/bruno/Documentos/DOC/R_Scripts_Numerical_Ecology/NEwR-2ed_code_data/NEwR2-Functions/sr.value.R")
source("/home/bruno/Documentos/DOC/R_Scripts_Numerical_Ecology/NEwR-2ed_code_data/NEwR2-Functions/quickMEM.R")
source("/home/bruno/Documentos/DOC/R_Scripts_Numerical_Ecology/NEwR-2ed_code_data/NEwR2-Functions/scalog.R")

load("~/Data100Lakes.RData")

# Data from flowDiv (100 Lakes)
mite=m1
mite.xy<-geoXY(dataall$Lat, dataall$Lon, unit=1000)
# Transform the data
mite.h <- mite
## dbMEM analysis of the 100 lakes data

# Is there a linear trend in the mite data?
anova(capscale(mite~., data.frame(mite.xy), dist="bray"))	# Result: significant trend for Bray-Curtis
## Step 1. Construct the matrix of dbMEM variables
mite.dbmem.tmp <- dbmem(mite.xy, silent = FALSE)
mite.dbmem <- as.data.frame(mite.dbmem.tmp)
# Truncation distance used above:
(thr <- give.thresh(dist(mite.xy)))

# Display and count the eigenvalues
attributes(mite.dbmem.tmp)$values
length(attributes(mite.dbmem.tmp)$values)
# Argument silent = FALSE allows the function to display 
# the truncation level.

## Step 2. Run the global dbMEM analysis on the *detrended*
##    Hellinger-transformed mite data
(mite.dbmem.rda <-capscale(mite.h~., mite.dbmem, dist="bray"))
anova(mite.dbmem.rda)


## Step 3. Since the R-square is significant, compute the adjusted
##    R2 and run a forward selection of the dbmem variables
(mite.R2a <- RsquareAdj(mite.dbmem.rda)$adj.r.squared)

# Descomentar para usar o ordistep
 mite.dbmem.ordistep <- ordistep(mite.dbmem.rda)

#(nb.sig.dbmem <- nrow(mite.dbmem.fwd))    # Number of signif. dbMEM
# Identity of the significant dbMEM in increasing order
#(dbmem.sign <- sort(mite.dbmem.fwd[ ,2]))
# Write the significant dbMEM to a new object
dbmem.red <- mite.dbmem[ ,c("MEM3", "MEM22","MEM1")]

## Step 4. New dbMEM analysis with 8 significant dbMEM variables
##    Adjusted R-square after forward selection: R2adj = 0.2418
(mite.dbmem.rda2 <-capscale(mite.h~., dbmem.red, dist="euclidian"))
#(mite.dbmem.rda2 <- rda(mite.h.det ~ ., data = dbmem.red))
(mite.fwd.R2a <- RsquareAdj(mite.dbmem.rda2)$adj.r.squared)
anova(mite.dbmem.rda2)
(axes.test <- anova(mite.dbmem.rda2, by = "axis"))
# Number of significant axes
(nb.ax <- length(which(axes.test[ ,ncol(axes.test)] <=  0.05)))

## Step 5. Plot the significant canonical axes
mite.rda2.axes <- 
  scores(mite.dbmem.rda2, 
         choices = c(1:nb.ax), 
         display = "lc", 
         scaling = 1)

par(mfrow = c(1,nb.ax))
for(i in 1:nb.ax){
  sr.value(mite.xy, mite.rda2.axes[ ,i], 
           sub = paste("RDA",i), 
           csub = 2)
}
# Interpreting the spatial variation: regression of the significant
# canonical axes on the environmental variables, with Shapiro-Wilk 
# normality tests of residuals
mite.env<-dataall[,c(9,13:37, 48:51)]
mite.rda2.axis1.env <- lm(mite.rda2.axes[ ,1] ~ ., data = mite.env)
shapiro.test(resid(mite.rda2.axis1.env))
summary(mite.rda2.axis1.env)

mite.rda2.axis2.env <- lm(mite.rda2.axes[ ,2] ~ ., data = mite.env)
shapiro.test(resid(mite.rda2.axis1.env))
summary(mite.rda2.axis1.env)

# Ptot e BGE (Bacterial Efficiency Growth?) influenciam os padrões espaciais
