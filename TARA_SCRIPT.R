
source("/home/brunomsw/Documentos/DOUTORADO/flowDiv/flowDiv/flowDiv/flowDiv_internals.R")
source("/home/brunomsw/Documentos/DOUTORADO/flowDiv/flowDiv/flowDiv/flowDiv_main.R")
library("vegan")
library("flowWorkspace")
library("ggplot2")
library("ggcyto")
library("plyr")
library("pangaear")
setwd("/home/brunomsw/Documentos/DOUTORADO/DATA_FCS/Tara/Results")

# Alpha  0.009419955
# Download Pangeaa
pangea.all=pg_data("10.1594/PANGAEA.842237")
pangea.data=pangea.all[[1]]$data
#pangea.data=data.frame(pangea.data)

############################## OM dataset
#OM_CompanionTable_MODIFIED <- read_delim("~/Documentos/DOUTORADO/DATA_FCS/Tara/TARA METADATA/OM.CompanionTable MODIFIED.csv", 
   #                                      +     "\t", escape_double = FALSE, locale = locale(decimal_mark = ","), 
    #                                     +     trim_ws = TRUE)

#Stations
#TARA_148b_MES_0.22-3 : search it before running

# First change
OM_CompanionTable_MODIFIED2=OM_CompanionTable_MODIFIED[!is.na(OM_CompanionTable_MODIFIED$OG.Shannon),]

# 
OM_regexpr=regexpr("[a-z]{4}_[0-9]{3}_[a-z]{3}", OM_CompanionTable_MODIFIED2$`Sample label [TARA_station#_environmental-feature_size-fraction]`, ignore.case = T)
OM_stations = regmatches(OM_CompanionTable_MODIFIED2$`Sample label [TARA_station#_environmental-feature_size-fraction]`,OM_regexpr, invert = F)

#Fractions

OM_fraction = regmatches(OM_CompanionTable_MODIFIED2$`Sample label [TARA_station#_environmental-feature_size-fraction]`,OM_regexpr, invert = T)
OM_fraction=factor(gsub("^_", "", do.call( rbind, OM_fraction)[,2]))

# Removing MIX 
mix=grep("MIX", OM_stations)
OM_CompanionTable_MODIFIED3=OM_CompanionTable_MODIFIED2[-mix,]
OM_stations=OM_stations[-mix]
OM_fraction=OM_fraction[-mix]

# New OM datasets
OM2=cbind(Station=OM_stations, Fractions=OM_fraction, OM_CompanionTable_MODIFIED3)

OM3=OM2[c(OM2$Fractions=="0.22-3"),]


############################## 
# flowDiv
tara.fd=flowDiv2("../Bruno_Tara_ALL_Bact10-3.wsp", "Bact", do.plot = F, use.beads = T,ialpha = "shannon", beads = "Beads", static = F, transform = F, autotrans=F, psize = 1)

id.flowdiv_raw=names(tara.fd$Alpha)
id.flowdiv_raw_regexpr=regexpr("[0-9]{1,3} [a-z]{1,4}", id.flowdiv_raw, ignore.case = T)

id.flowdiv_clean=regmatches(id.flowdiv_raw,id.flowdiv_raw_regexpr, invert = F)
# Separates numbers and locals from flowDiv (to fit Pangea)
## Numbers
fd.numbers=regexpr("[0-9]{1,3}", id.flowdiv_clean, ignore.case = T)
fd.numbers1 = regmatches(id.flowdiv_clean,fd.numbers)

# Adding "zeros" # Check samples  87, 117, 118 before running
# One char
ones=nchar(fd.numbers1)==1
fd.numbers1[ones]=paste("00",sep="", fd.numbers1[ones])
# Two chars
twos=ones=nchar(fd.numbers1)==2
fd.numbers1[twos]=paste("0",sep="", fd.numbers1[twos])

#Fitting Pangea tara1
tara1=paste("TARA_", fd.numbers1, sep="") # Use tara1 for Panagea dataset
tara1=data.frame(Station=tara1)
## Locals
fd.letters=regexpr("[a-z]{1,4}", id.flowdiv_clean, ignore.case = T)
fd.letters1= regmatches(id.flowdiv_clean,fd.letters)
# Fiting OM dataset tara2
fd.letters1[grep("SU", fd.letters1, ignore.case = T)]="SRF"
fd.letters1[grep("ME", fd.letters1, ignore.case = T)]="MES"
fd.letters1[grep("DCM", fd.letters1, ignore.case = T)]="DCM"
tara2=paste(as.matrix(tara1),"_", fd.letters1, sep="") 
tara2=data.frame(tara2)
names(tara2)="Station"



ids=tara2$Station%in%OM3$Station

#  Enviromental matrices

#tara.env_PANGEA=join(tara1, pangea.data, by="Station")

OM4=OM3[ids,]
tara3=data.frame("Station"=tara2[ids,])
tara.env_OMG=join(tara3, OM3, by="Station")

cor1=corr.test(cbind(Index=tara.fd$Alpha[ids], tara.env_OMG)[c(1,29:42)])


# Teste 01

m1=tara.fd$Matrices
rownames(m1)

ord <- metaMDS(m1, try=100)
(fit <- envfit(ord, tara.env_PANGEA[,c(5,6,12,13,14)], perm = 999, na.rm = T))
scores(fit, "vectors")

plot(ord, type="n")
points(ord$points)
plot(fit, p.max = 0.05, col = "red")

# Teste 02
alphas=cbind(tara2, Alpha=tara.fd$Alpha)
pielous=cbind(tara2, Pielou=tara.fd$Pielou)

join1=join(OM2,alphas, by="Station")
join2=join(join1,pielous, by="Station")


#Teste 03
OM3=OM2[c(OM2$Fractions=="0.22-3"),]
t1=join(alphas, OM3, by="Station")
m1=tara.fd$Matrices

"shannon"

best.index<-function(w,ind="NULL"){
  if(ind=="NULL") m1.shannon=apply(m1, 1, function(x) diversity(x,index = w))
  if(ind=="rich") m1.shannon=apply(m1, 1, function(x) specnumber(x, MARGIN = 1))
  if(ind=="pielou")m1.shannon<-apply(m1, 1, function(x)diversity(x)/log(vegan::specnumber(x)))
  m2.shannon=cbind(tara2, m1.shannon)
  
  t1=join( m2.shannon, OM3, by="Station")
  t2=t1[!is.na(t1$Fractions),]
  corre=corr.test(t2[,c(2,c(8:42))])
   m1=cbind("p"=data.frame(corre$p[,1]), "r"=data.frame(corre$r[,1]))
   m2=m1[m1[,1]<=0.05,]
  return(m2)
}


best.index("shannon") # *

best.index("invsimpson")
best.index("simpson")

best.index("shannon", ind="rich")
best.index("shannon", ind="pielou")

#correlation(package psych)
corre=corr.test(t2[,8:43])
corre$p
