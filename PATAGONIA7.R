## Required Packages
library("betapart")
library("reshape")
library("vegan")
library("flowWorkspace")
library("ggplot2")
library("ggcyto")
library("plyr")
library("car")
library("corrplot")
library("pangaear")
library("ggbiplot")
library("psych")
library("cowplot")
library("RColorBrewer")
library("dunn.test")
library("seriation")
library("ggbiplot")
library("gvlma")
library("RVAideMemoire")
library("gplots")
library("gridExtra")
## Sourcing flowDiv
source("/home/bruno/Documentos/DOC/flowDiv_SUBMISSION/Packages_and_Codes/flowDiv2/R/flowDiv_internalsCORRECTED.R")
source("/home/bruno/Documentos/DOC/flowDiv_SUBMISSION/Packages_and_Codes/flowDiv2/R/flowDiv_main1.R")
## Sourcing data
source("/home/bruno/Documentos/DOC/MAIN/Patagonia_OK/Scripts_Data/DATA_IMPORT.R")
## Workings directly
setwd("/home/bruno/Documentos/DOC/MAIN/Patagonia_OK/Results_AFTER_PLOS")

# Removing Huapi and P.10
env<-env[-4,]
dgge<-dgge[-c(4,10),]
# Correcting ID
#env
env$`Water body`[24]="SL. Negra"
env$`Water body`[31]="L. Yehuín"
env$`Water body`[21]="SL. de los Cisnes"
#DATA
DATA$Lago[13]="P.15"
DATA$Lago[21]="L. Acigami"
DATA$Lago[28]="SL. Verde"
## Summary trophic status
table(DATA$`TROPHIC STATE`)
### Merging into meso-eutrofic
DATA$`TROPHIC STATE`=revalue(DATA$`TROPHIC STATE`, c("EUT"="MESO_EUT", "MESO"="MESO_EUT", "OLIGO"="OLIGO"))
# dgge
dgge$X1[22]="SL. Victoria"
dgge$X1[27]="L. Yehuín"
dgge$X1[26]="SL. de los Cisnes"
# Correcting blank spaces
env$`Water body`=gsub("\\.","\\. ",gsub("\\. ", "\\.", env$`Water body`),env$`Water body`)
DATA$Lago=gsub("\\.","\\. ",gsub("\\. ", "\\.", DATA$Lago),DATA$Lago)
dgge$X1=gsub("\\.","\\. ",gsub("\\. ", "\\.", dgge$X1),dgge$X1)
# Visual checking 1
cbind(DATA$Lago, env$`Water body`)
# Matching sites
match(DATA$Lago, env$`Water body`)
# New ENV
env2=env[match(DATA$Lago, env$`Water body`),]
# Visual checking 2
cbind(DATA$Lago, env2$`Water body`)

# flowDiv
# 3 dimensions, 75 bins per dimesnion
#
#pata.fd=flowDiv("/home/bruno/Documentos/DOC/MAIN/Patagonia_OK/METADATA/flowDiv patagonia/Romina_Bruno_1_10.wsp", "Bact tot", do.plot = F, static = F, dilutions=DATA$`x (flowDiv)`, transform = T, autotrans=T,psize = 1.08)
# 2 dimensions
#pata.fd=flowDiv("/home/bruno/Documentos/DOC/MAIN/Patagonia_OK/METADATA/flowDiv patagonia/Romina_Bruno_1_10.wsp", "Bact tot", do.plot = T, static = F, dilutions=DATA$`x (flowDiv)`, transform = T, autotrans=T,psize = 1.08)


# Getting count

#wksp=opc("/home/bruno/Documentos/DOC/MAIN/Patagonia_OK/METADATA/flowDiv patagonia/Romina_Bruno_1_10.wsp")
#samples=nn(wksp,nod = "Bact tot", use.beads = F, nod2="NULL")
counts=unlist(lapply(wksp[[1]], function(x)getTotal(x, "Bact tot")))
counts.corrected=counts*DATA$`x (flowDiv)`



# ID's flowDiv
id.flowdiv_raw=names(pata.fd$Alpha)
id.flowdiv_raw_regexpr=regexpr("[0-9]{2} dil [0-9]{1}_[0-9]{1}", id.flowdiv_raw, ignore.case = T)
id.flowdiv_clean=regmatches(id.flowdiv_raw,id.flowdiv_raw_regexpr, invert = F)

# ID's DATA

id.data_raw=DATA$file
id.data_raw_regexpr=regexpr("[0-9]{2} dil [0-9]{1}_[0-9]{1}", id.data_raw, ignore.case = T)
id.data_clean=regmatches(id.data_raw,id.data_raw_regexpr, invert = F)
DATA2=cbind("ID"=id.data_clean, DATA)

# Visual checking 3
cbind(DATA$Lago, env2$`Water body`, DATA2$Lago)
(DATA$Lago==env2$`Water body`)&&(env2$`Water body`== DATA2$Lago)
# Changing names
rownames(pata.fd$Matrices)<-as.character(DATA2$Lago)

# Matriz renomeada
m1=pata.fd$Matrices

#Definindo os indices
rich=apply(m1, 1, function(x) specnumber(x, MARGIN = 1))
pielou=apply(m1, 1, function(x)diversity(x)/log(vegan::specnumber(x)))
shannon=apply(m1, 1, function(x) diversity(x,index = "shannon"))
simpson=apply(m1, 1, function(x) diversity(x,index = "simpson"))
invsimpson=apply(m1, 1, function(x) diversity(x,index = "invsimpson"))

# New data set

DATA3=cbind(DATA2, "FC_Richness"=rich, "FC_Pielou"=pielou, "FC_Shannon"=shannon, "FC_Simpson"=simpson, "FC_invSimpson"=invsimpson)
# Adjusting labesl
label1=gsub("MESO_EUT", "MESO-EUTROPHIC", DATA3$`TROPHIC STATE`)
label1=gsub("OLIGO", "OLIGOTROPHIC", label1)
# Plot
mylabels=factor(label1)
# Wilcox
DATA4=DATA3
colnames(DATA4)[7]<-"TROPHIC_STATE"
# Salvo daqui para cima!
#wilcox.test(FC_Richness~TROPHIC_STATE, DATA4)
wilcox.test(log10(FC_Richness)~TROPHIC_STATE, DATA4) # significant
wilcox.test(FC_Pielou~TROPHIC_STATE, DATA4)
wilcox.test(FC_Shannon~TROPHIC_STATE, DATA4)
wilcox.test(FC_Simpson~TROPHIC_STATE, DATA4)
wilcox.test(FC_invSimpson~TROPHIC_STATE, DATA4)
# Boxplots

DATA3df=data.frame(DATA3)
DATA3df$Log_Richness=log10(DATA3$FC_Richness)

mybox<-function(x, lab){
  g0=ggplot(aes_string(y=x, x="TROPHIC.STATE"), data=DATA3df)+geom_boxplot(aes_string(fill=mylabels, alpha=.5))+ theme(legend.position="none")+ scale_x_discrete(name="")+geom_point(size=3)+
    theme(text = element_text(size=20),
          axis.line=element_blank(),
          axis.ticks=element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_text(size=20))+
    labs(y = lab)
  return(g0)
}


b1=mybox("FC_Shannon", "Shannon Index")
b2=mybox("FC_Pielou", "Pielou Index")
#mybox("FC_Richness", "Richness")
b3=mybox("Log_Richness", expression('Log'[10]*' Richness'))
# Density
b4=ggplot(DATA3df, aes(x=Log_Richness))+
  geom_density(aes(group=TROPHIC.STATE, colour=TROPHIC.STATE, fill=TROPHIC.STATE), alpha=0.3)+
  theme(text = element_text(size=20),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        axis.text.y=element_text(size=20),
        legend.position="none")+
  labs(y = "Density", x=expression('Log'[10]*' Richness'))


grid.arrange(b1,b2,b3,b4)


# Correlations

cor.mtest <- function(mat, conf.level = 0.95){
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      tmp <- cor.test(mat[,i], mat[,j], conf.level = conf.level, method = "spearman", exact = F)
      p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
      
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}

M4=cbind(FC_Shannon=DATA3$FC_Shannon,FC_Pielou=DATA3$FC_Pielou, FC_Richness=DATA3$FC_Richness,env2[,-c(1:3)])


corr.test(M4, method = "spearman")
res1 <- cor.mtest(M4, 0.95)
# Adjusting colnames
M4.1=M4
colnames(M4.1)<-c("Shannon","Pielou","Richness","Latitude","Longitude","Altitude","Area","Temperature","pH","Conductivity","DO","DIN","Kd","Chl a","Phosphate","DOC")
corrplot(cor(M4.1, method = "spearman", use="complete.obs"), p.mat = res1[[1]], sig.level=.05, type = "lower", pch.col="gray20", method = "square", cl.cex=.8, number.cex = .8,
         col=c("red3", "blue3"), tl.col = "black", pch.cex = 2.5)


# Scatterplots

M5=M4
colnames(M5)
colnames(M5)[10]="pH"
colnames(M5)[9]="Temperature"
colnames(M5)[15]="Chla"
M5$Chla = log10(M5$Chla)

mycor<-function(x,y){
  g0=ggplot(aes_string(x, "DATA3$FC_Shannon"), data=M5)+geom_text(aes(label=`Water body`))+geom_smooth(method='lm')+
    xlab(y)+ylab("Shannon Index")
  return(g0)
}

# Richness versus Shannon
ggplot(aes_string("FC_Richness", "FC_Shannon"), data=DATA3df)+geom_text(aes(label=Lago))+geom_smooth(method='lm')+
  xlab("FC_Richness")+ylab("Shannon Index")

ggplot(aes(log10(FC_Richness), FC_Shannon), data=DATA3df)+geom_point()+geom_smooth(method='lm')+
  xlab("FC_Richness")+ylab("Shannon Index")


# Other plots
mycor("pH", "pH")
mycor("Temperature", "Temp.(°C)")
mycor("Chla", "Log10 Chl a (μg L−1)")


# PCA

M5=M4
colnames(M5)[1:3]<-c("Shannon",  "Pielou",   "Richness")
wine.pca <- prcomp(M5[,c(1:3)], scale. = TRUE)
pts<-data.frame(wine.pca$x)
vects<-data.frame(wine.pca$rotation)

summary(wine.pca)
ggplot(aes(PC1, PC2, colour=mylabels), data=pts)+
  geom_point(size=5)+
  stat_ellipse(level=0.95, geom = "polygon", alpha = .2, aes(fill = mylabels))+
  geom_segment(data=vects, mapping=aes(x=0, y=0, xend=PC1, yend=PC2), arrow=arrow(length = unit(0.2, "cm")), size=.8, color="black", linetype=1)+
  geom_text(data=vects, aes(x=PC1*1.2, y=PC2*1.2, label=rownames(vects)), size=8, vjust=0,hjust=c(0,0,0), colour="black")+
  theme_bw()+
  theme(panel.border = element_blank(),
        text = element_text(size=23),
        legend.title=element_blank(),
        legend.background = element_rect(colour="white", size=.5, linetype=3),
        legend.position = "bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )+
  labs(x = "PC1 (65.59%)", y = "PC2(32.85%)")


meta2=metaMDS(vegdist(m1), autotransform = T, try = 200)

plot1=function(x){
  
  plot(x, type="n")
  text(x$points, labels=DATA3$Lago, col=as.numeric(DATA3$`TROPHIC STATE`), pch=17)
  
}

plot1(meta2)


meta2=metaMDS(m1, autotransform = T, try = 200)
fit=envfit(meta2, M5, na.rm = T)
fit
variables<-rownames(fit$vectors$arrows)[fit$vectors$pvals<=0.05]


plot(meta2, display="sites")
plot(fit, p.max=0.05)

for (i in 1:12){
  ordisurf(meta2~M5[,variables[i]], main=variables[i])
  arrows(0,0, fit$vectors$arrows[variables[i],1],fit$vectors$arrows[variables[i],2])
} 


fit3<-fit$vectors$arrows[fit$vectors$pvals<=0.05,]
(fit3=data.frame(fit3))
rownames(fit3)<-c("Shannon","Pielou","Richness","Lat.","Lon.","Area","Temp.","pH","Cond.","Kd","Chla","DOC")

fit3=fit3[c("Lat.","Lon.","Area","Temp.","pH","Cond.","Kd","Chla","DOC"),]
# # Visual adjust
# fit3=fit2[c(5:10,15,17),]*c(rep(1.1, 6), 1.2, 1.2)
# fit3[8,1]=fit3[8,1]+.6

ggplot(aes(NMDS1, NMDS2, colour=mylabels), data=data.frame(scores(meta2)))+
  geom_point(size=2)+
  stat_ellipse(level=0.95, geom = "polygon", alpha = .2, aes(fill = mylabels))+
  geom_segment(data=fit3, mapping=aes(x=0, y=0, xend=NMDS1, yend=NMDS2), arrow=arrow(length = unit(0.3, "cm")), size=.6, color="black", linetype=1)+
  geom_segment(data=fit3[c(7,9),], mapping=aes(x=0, y=0, xend=NMDS1*.95, yend=NMDS2*.95), arrow=arrow(length = unit(0.3, "cm")), size=.6, color="black", linetype=1)+
  geom_text(data=fit3, aes(x=NMDS1, y=NMDS2), label=rownames(fit3), size=5, vjust=c(.5,-.2,0,1,1,1,2,1,2),hjust=c(0,1,0,1,1,-.5,0,1,.5), colour="black")+
  theme_bw()+
  theme(panel.border = element_blank(),
        text = element_text(size=15),
        legend.title=element_blank(),
        legend.background = element_rect(colour="white", size=.5, linetype=3),
        legend.position = "bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )


# Adonis
adonis(vegdist(m1)~DATA3$`TROPHIC STATE`+env2$`Long. (°W)`+env2$`Lat. (°S)`)

# Beta disper # Not significant!

anova(betadisper(vegdist(m1), DATA3$`TROPHIC STATE`))

# Pairwise
pairwise.perm.manova(vegdist(m1),DATA3$`TROPHIC STATE`, p.method = "bonferroni")

# New variables
m2=vegdist(m1)
m3=as.matrix(m2)

meso= DATA3$`TROPHIC STATE`=="MESO_EUT"
oligo= DATA3$`TROPHIC STATE`=="OLIGO"


m=m3[meso,meso]
o=m3[oligo,oligo]


o0=o[upper.tri(o)]
m0=m[upper.tri(m)]

all=stack(list("o"=c(o), "m"=c(m)))
all0=stack(list("o"=c(o0), "m"=c(m0)))


boxplot(values~ind, all, ylim=c(0,1))
boxplot(values~ind, all0, ylim=c(0.3,1))

# Bray-Curtis Distance
ggplot(aes(y=values, x=ind), data=all0)+geom_boxplot(aes(fill=ind, alpha=.5))+
  theme(legend.position="none",
        text = element_text(size=25),
        axis.text=element_text(size=25),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank())+
  labs(y ="Bray-Curtis Distance")

# Variance

a=betadisper(vegdist(m1), DATA3$`TROPHIC STATE`)
df=data.frame(a$distances, a$group)

ggplot(aes(y=a.distances, x=a.group), data=df)+geom_boxplot(aes(fill=a.group, alpha=.5))+
  theme(legend.position="none",
        text = element_text(size=25),
        axis.text=element_text(size=25),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank())+
  labs(y ="Distance to Centroid")


# Turnover x Nestdness
oligo=apply(m1[oligo,], 2, sum)/nrow(m1[oligo,])
meso=apply(m1[meso,], 2, sum)/nrow(m1[meso,])
all=rbind("oligo"=oligo, "meso"=meso)
mall=vegdist(all)
braym=bray.part(all)
nest=braym$bray.bal/braym$bray
turn=1-nest

df1=data.frame(cbind("contrast"=c("MxO"), "TURNOVER"=c(turn), "NESTEDNESS"=c(nest)))
df2=melt(df1, id.var="contrast")
df2$value=as.numeric(as.character(df2$value))

pie<-ggplot(df2, aes(x = "", y = value, fill = variable))+
  geom_bar(stat = "identity", alpha=0.9)+
  coord_polar("y", start=0)
  
pie+scale_y_continuous(labels = scales::percent)+
  theme(text = element_text(size=13),
        axis.line = element_blank(),
        axis.text=element_text(size=10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        legend.position = "top")+
scale_fill_manual(values=c("mediumvioletred", "cyan4"))
 # labs(title = "Bray-Curtis Distance (%)")

