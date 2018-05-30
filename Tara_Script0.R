source("/home/brunomsw/Documentos/DOUTORADO/flowDiv/flowDiv/flowDiv/flowDiv_internals.R")
source("/home/brunomsw/Documentos/DOUTORADO/flowDiv/flowDiv/flowDiv/flowDiv_main.R")
library("vegan")
library("flowWorkspace")
library("ggplot2")
library("ggcyto")
library("plyr")
setwd("/home/brunomsw/Documentos/DOUTORADO/DATA_FCS/Tara/Results")




# flowDiv
t1=flowDiv2("../Bruno_Tara_ALL_Bact10-3.wsp", "Bact", do.plot = T, use.beads = T, beads = "Beads", static = F, transform = F, autotrans=F, psize = 1)

# Regex
id5=names(t1$Alpha)
id5=gsub(" ","", id5)
a=regexpr("[a-z]{2,4}[0-9]{1,3}[a-z]{1,4}", id5, ignore.case = T)
b=regmatches(id5,a)
c=gsub("^[a-z]{2,4}","", b, ignore.case = T)

j=grep("^[0-9]{1}[a-z]", c, ignore.case = T)
c[j]=paste("0",sep="", c[j])
l=grep("^[0-9]{2}[a-z]", c, ignore.case = T)
c[l]=paste("0",sep="", c[l])


c=gsub("su[a-z]{1,2}", "_SRF", c, ignore.case = T)
c=gsub("me[a-z]{1,2}", "_MES", c, ignore.case = T)
c=gsub("dc[a-z]{1,2}", "_DCM", c, ignore.case = T)



a2=regexpr("[0-9]{3}_[A-Z]{3}", labs[,1], ignore.case = T)
b2=regmatches(labs[,1],a2)
c2=gsub("^[a-z]{2,4}","", b2, ignore.case = T)




m1=metaMDS(t1$Matrices, try = 100)
plot(m1, type="n")
text(m1$points, labels = c, cex=.5)
km=kmeans(m1$points, 2)
c[km$cluster==1]
