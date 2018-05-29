#########################################################################
### Code for statistical analysis of flowDiv
### Bruno Mattos S. Wanderley, Universidade Federal do RN, Brazil
#########################################################################
##################### Required Packages #################################

####################### Sourcing codes ##################################
source("/home/bruno/Documentos/DOC/flowDiv_SUBMISSION/Packages_and_Codes/flowDiv2/R/flowDiv_internalsCORRECTED.R")
source("/home/bruno/Documentos/DOC/flowDiv_SUBMISSION/Packages_and_Codes/flowDiv2/R/flowDiv_main1.R")
####################### Sourcing data ###################################
source("/home/bruno/Documentos/DOC/MAIN/Patagonia_OK/Scripts_Data/DATA_IMPORT.R")
## Sourcing data
# DATA_IMPORT.R contains: dgge, env DATA and Cybar_counts
# dgge = DGGE band's information
# env = enviromental variables
# DATA = ID's and dilutions for samples
# Cybar_counts = gates counts used in FlowCybar pipeline
##################### Workings directory ################################
setwd("/home/bruno/Documentos/DOC/MAIN/Patagonia_OK/Results_AFTER_PLOS")
##################### Data preparation 1 ################################
library("flowWorkspace")
library("plyr")
library("vegan")
## Removing Huapi and P.10 (Reason: Low quality cytograms)
env<-env[-4,]
dgge<-dgge[-c(4,10),]
## Correcting ID's (Data came from diferent files. Threre are minor misspellings that can impair correct matching)
# env
env$`Water body`[24]="SL. Negra"
env$`Water body`[31]="L. Yehuín"
env$`Water body`[21]="SL. de los Cisnes"
# DATA
DATA$Lago[13]="P.15"
DATA$Lago[21]="L. Acigami"
DATA$Lago[28]="SL. Verde"
# dgge
dgge$X1[22]="SL. Victoria"
dgge$X1[27]="L. Yehuín"
dgge$X1[26]="SL. de los Cisnes"
# Correcting blank spaces
env$`Water body`=gsub("\\.","\\. ",gsub("\\. ", "\\.", env$`Water body`),env$`Water body`)
DATA$Lago=gsub("\\.","\\. ",gsub("\\. ", "\\.", DATA$Lago),DATA$Lago)
dgge$X1=gsub("\\.","\\. ",gsub("\\. ", "\\.", dgge$X1),dgge$X1)
# Checking matching
sum(!(DATA$Lago==env$`Water body`))!=0
# Correctly matching sites
match1<-match(DATA$Lago, env$`Water body`)
# New env (used only to not override previous dataset)
env2=env[match1,]
# Confirming matching
sum(!(DATA$Lago==env2$`Water body`))!=0

## Checking balance of data set
# Trophic status counts
table(DATA$`TROPHIC STATE`)
# Merging into meso-eutrofic (Reason: threre is a big imbalance between factors)
DATA$`TROPHIC STATE`=revalue(DATA$`TROPHIC STATE`, c("EUT"="MESO_EUT", "MESO"="MESO_EUT", "OLIGO"="OLIGO"))

##################### Data from flowDiv #################################
# Run with the folloing parameters: 
# Channels: FL1-H, FL3-H and SSC-C (used previously for sequential gating on FlowJo)
# Number os bins: 75 (as suggested by flowDiv)
# Number of clusters: 5 
pata.fd=flowDiv("/home/bruno/Documentos/DOC/MAIN/Patagonia_OK/METADATA/flowDiv patagonia/Romina_Bruno_1_10.wsp", "Bact tot", do.plot = F, static = F, dilutions=DATA$`x (flowDiv)`, transform = T, autotrans=T,psize = 1.08)
2 3 5
75
#5
## Getting flowSet only (for comparisons)
wksp=opc("/home/bruno/Documentos/DOC/MAIN/Patagonia_OK/METADATA/flowDiv patagonia/Romina_Bruno_1_10.wsp")
samples=nn(wksp,nod = "Bact tot", use.beads = F, nod2="NULL")
##################### Data preparation 2 ################################
## Reformating ID's from flowDiv and DATA to check matching between sites

# ID's flowDiv
id.flowdiv_raw=names(pata.fd$Alpha)
id.flowdiv_raw_regexpr=regexpr("[0-9]{2} dil [0-9]{1}_[0-9]{1}", id.flowdiv_raw, ignore.case = T)
id.flowdiv_clean=regmatches(id.flowdiv_raw,id.flowdiv_raw_regexpr, invert = F)

# ID's DATA

id.data_raw=DATA$file
id.data_raw_regexpr=regexpr("[0-9]{2} dil [0-9]{1}_[0-9]{1}", id.data_raw, ignore.case = T)
id.data_clean=regmatches(id.data_raw,id.data_raw_regexpr, invert = F)

# Checking matching
sum(!(id.data_clean==id.flowdiv_clean))==0
# New data set with ID's column reformated and checked
DATA2=cbind("ID"=id.data_clean, DATA)
# Just checking
(DATA$Lago==env2$`Water body`)&&(env2$`Water body`== DATA2$Lago)
# Changing the row names of raw matrix from flowDiv
rownames(pata.fd$Matrices)<-as.character(DATA2$Lago)

# New matrix m1 to store renamed raw matrix from flowDiv
m1=pata.fd$Matrices

# Defining diversity index for m1. 
# Note: One can get those indices from flowDiv directly or by this way.
rich=apply(m1, 1, function(x) specnumber(x, MARGIN = 1))
pielou=apply(m1, 1, function(x)diversity(x)/log(vegan::specnumber(x)))
shannon=apply(m1, 1, function(x) diversity(x,index = "shannon"))

# New data set to store the new indices created previously
DATA3=cbind(DATA2, "FC_Richness"=rich, "FC_Pielou"=pielou, "FC_Shannon"=shannon)
# Renaming labels for trophic status
label1=gsub("MESO_EUT", "MESO-EUTROPHIC", DATA3$`TROPHIC STATE`)
label1=gsub("OLIGO", "OLIGOTROPHIC", label1)
# Changing classes to factor
mylabels=factor(label1)

######################### Statistics ####################################
library("gvlma")
# New data set DATA4 (again, just to not overwrite)
DATA4=DATA3
colnames(DATA4)[7]<-"TROPHIC_STATE"

# Checking adequacy of models for ANOVA
gvlma(lm(DATA4$FC_Shannon~DATA4$TROPHIC_STATE)) #OK
gvlma(lm(DATA4$FC_Pielou~DATA4$TROPHIC_STATE)) #OK
gvlma(lm(log10(DATA4$FC_Richness)~DATA4$TROPHIC_STATE)) # OK

# Running parametric ANOVA
summary(aov(lm(DATA4$FC_Shannon~DATA4$TROPHIC_STATE))) # Not significant
summary(aov(lm(DATA4$FC_Pielou~DATA4$TROPHIC_STATE))) # Not significant
summary(aov(lm(log10(DATA4$FC_Richness)~DATA4$TROPHIC_STATE))) # Significant


# Wilcoxon test (not used as assumptions for paramteric ANOVA are not violeted)
wilcox.test(FC_Shannon~TROPHIC_STATE, DATA4) # Not significant
wilcox.test(FC_Pielou~TROPHIC_STATE, DATA4) # Not significant
wilcox.test(log10(FC_Richness)~TROPHIC_STATE, DATA4) # Significant

############################ Plots ######################################
library("corrplot")
library("ggcorrplot")
library("ggplot2")
library("gridExtra")
library("psych")
### New functions for plotting
## For boxplots
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
## For correlations
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

# Converting DATA3 to dataframe
DATA3df=data.frame(DATA3)
DATA3df$Log_Richness=log10(DATA3$FC_Richness)

## Plots
# Boxplots
b1=mybox("FC_Shannon", "Shannon Index")
b2=mybox("FC_Pielou", "Pielou Index")
b3=mybox("Log_Richness", expression('Log'[10]*' Richness'))
# Density plot
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

# Plotting all in a 2x2 grid
grid.arrange(b1,b2,b3,b4)


# Correlations
# M4 to concatenate data
M4=cbind(FC_Shannon=DATA3$FC_Shannon,FC_Pielou=DATA3$FC_Pielou, FC_Richness=DATA3$FC_Richness,env2[,-c(1:3)])
# Correlations based on Spearman method
corr.test(M4, method = "spearman")
res1 <- cor.mtest(M4, 0.95)
# Adjusting colnames
colnames(M4)<-c("Shannon","Pielou","Richness","Latitude","Longitude","Altitude","Area","Temperature","pH","Conductivity","DO","DIN","Kd","Chla","Phosphate","DOC")
# Log10 for Chl a
M4$Chla = log10(M4$Chla)
# Correlation plot
corrplot(cor(M4, method = "spearman", use="complete.obs"), p.mat = res1[[1]], sig.level=.05, type = "upper", pch.col="gray20", method = "square", cl.cex=.8, number.cex = .8,
         col=c("red3", "blue3"), tl.col = "black", pch.cex = 2.5)
#################
# Correlations
res2<-cor_pmat(M4, method = "spearman", use="complete.obs")
# Visual checking
round(res2, 2)
ggcorrplot(cor(M4, method = "spearman", use="complete.obs"), 
           p.mat = res2, sig.level=.05, type = "upper", 
           legend.title = "")
#################
# Scatterplots
library("grid")
g1=ggplot(aes(log10(Richness), Shannon), data=M4)+geom_point()+geom_smooth(method='lm')+
  labs(x=expression(Log[10]~Richness), y="")

g2=ggplot(aes(pH, Shannon), data=M4)+geom_point()+geom_smooth(method='lm')+
  labs(y="")

g3=ggplot(aes(log10(Kd), Shannon), data=M4)+geom_point()+geom_smooth(method='lm')+
  labs(x=expression(Log[10]~Kd), y="")

g4=ggplot(aes(log10(DOC), Shannon), data=M4)+geom_point()+geom_smooth(method='lm')+
  labs(x=expression(Log[10]~DOC), y="")


grid.arrange(arrangeGrob(g1 + theme(legend.position="none"), 
                         g2 + theme(legend.position="none"),
                         g3 + theme(legend.position="none"),
                         g4 + theme(legend.position="none"), 
                         nrow = 2,
                         top = textGrob("", vjust = 1, gp = gpar(fontface = "bold", cex = 1.5)),
                         left = textGrob("Shannon Index", rot = 90, vjust = 1, gp = gpar(fontface = "bold", cex = 1.5))))


################## PCA
wine.pca <- prcomp(M4[,c(1:3)], scale. = TRUE)
pts<-data.frame(wine.pca$x)
vects<-data.frame(wine.pca$rotation)
# Summary of PCA
summary(wine.pca)
# Plotting PCA
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
  )+ labs(x = "PC1 (65.59%)", y = "PC2(32.85%)")
################## NMDS

meta2=metaMDS(m1, autotransform = T, try = 200)
fit=envfit(meta2, M4, na.rm = T)
# Storing variables with pvals <= 0.05
variables<-rownames(fit$vectors$arrows)[fit$vectors$pvals<=0.05]

# Default plotting 
plot(meta2, display="sites")
plot(fit, p.max=0.05)
# Visual checking of linearity of variables
for (i in 1:length(variables)){
  ordisurf(meta2~M4[,variables[i]], main=variables[i])
  arrows(0,0, fit$vectors$arrows[variables[i],1],fit$vectors$arrows[variables[i],2])
} 

# Plot withh ggplot
fit3<-fit$vectors$arrows[fit$vectors$pvals<=0.05,]
# Removing diversity indices from plot
(fit3=data.frame(fit3[-c(1:3),]))
# Plot
ggplot(aes(NMDS1, NMDS2, colour=mylabels), data=data.frame(scores(meta2)))+
  geom_point(size=2)+
  stat_ellipse(level=0.95, geom = "polygon", alpha = .2, aes(fill = mylabels))+
  geom_segment(data=fit3, mapping=aes(x=0, y=0, xend=NMDS1, yend=NMDS2), arrow=arrow(length = unit(0.3, "cm")), size=.6, color="black", linetype=1)+
  geom_segment(data=fit3[9,], mapping=aes(x=0, y=0, xend=NMDS1*.95, yend=NMDS2*.95), arrow=arrow(length = unit(0.3, "cm")), size=.6, color="black", linetype=1)+
  geom_text(data=fit3, aes(x=NMDS1, y=NMDS2), label=rownames(fit3), size=5, vjust=c(.5,-.2,0,1,1,1,2,1,2.2),hjust=c(0,1,0,1,1,.5,0,1,.5), colour="black")+
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

# PERMANOVA accounting for spatial variation
adonis(vegdist(m1)~DATA3$`TROPHIC STATE`+env2$`Long. (°W)`+env2$`Lat. (°S)`)
# Homogeinity of groups assumptions test # Not significant
anova(betadisper(vegdist(m1), DATA3$`TROPHIC STATE`))

################## Nestdeness and Turnover
library("betapart")
library("reshape")
# Grouping
meso= DATA3$`TROPHIC STATE`=="MESO_EUT"
oligo= DATA3$`TROPHIC STATE`=="OLIGO"

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
# Raw plot
pie<-ggplot(df2, aes(x = "", y = value, fill = variable))+
  geom_bar(stat = "identity", alpha=0.9)+
  coord_polar("y", start=0)
# Complete plot
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

#########################################################################
############################# CHIC #####################################
# Link to original code:
# http://www.ufz.de/index.php?en=38441
#########################################################################
# Reading files
#reads overlapt results estimated by ImageJ
overlaps<-read.delim("/home/bruno/Documentos/DOC/KOCH_PIPELINES/CHIC/Plots/Results_overlaps.txt")     
#reads XOR results estimated by ImageJ
xor<-read.delim("/home/bruno/Documentos/DOC/KOCH_PIPELINES/CHIC/Plots/Results_xor.txt")      
#removes "overlap_" from file labels
overlaps$Label<-gsub("overlap_","",overlaps$Label)    
#removes ".bmp" from file labels
overlaps$Label<-gsub(".bmp","",overlaps$Label)     
#crates a file with proper label names
label.split<-strsplit(as.character(overlaps$Label),"_&_")  
#creates an empty vector file
nam<-NULL
#creates a file ("nam") with all combinations of label names in it 
for (i in 1:length(label.split)) nam<-rbind(nam, cbind(label.split[[i]][1], label.split[[i]][2])) 

dats<-data.frame(Label=overlaps$Label, Label.1=nam[,1], Label.2=nam[,2], IntDen=xor$IntDen, Area=overlaps$Area, IntDen.Area=xor$IntDen/overlaps$Area/100)      
#number of unique data labels is estimated
size<-length(unique(dats$Label.1))+1
#empty matrix ('mat') for disimilariy matrix is created
mat<-matrix(nrow=size, ncol=size)      
#counter J created and set to 1
j=1          
#dissimilarity according to formula 1 in MS for each combination is calculated and stored in matrix 'mat'
for (i in unique(dats$Label.1)) {      
  mat[,j]<-c(rep(NA,size-length(dats$IntDen.Area[dats$Label.1==i])),dats$IntDen.Area[dats$Label.1==i])
  j<-j+1
}
#diagonals are filled up with '0' values
diag(mat)<-0           
#labels names are written into matrix
colnames(mat)<-rownames(mat)<-union(unique(dats$Label.1), unique(dats$Label.2)) 
#matrix mat is converted to dissimilarity matrix ('mat.dist') to be used with metaMDS
mat.chic<-as.dist(mat)        
#Raw matrix
mat.chic2<-as.matrix(mat.chic)

# Renamed matrix
colnames(mat.chic2)<-rownames(mat.chic2)<-DATA3$Lago

#########################################################################
############################# flowCyBar #################################
# Link to original code:
# http://www.bioconductor.org/packages/devel/bioc/vignettes/flowCyBar/inst/doc/flowCyBar-manual.pdf
#########################################################################
library("flowCyBar")
cybar=flowCyBar::normalize(Cybar_counts[,-1],digits=2)

# Correcting problems for vegan (class of cybar is "AsIs")
cybar$B1=as.numeric(cybar$B1)
cybar$B2=as.numeric(cybar$B2)
cybar$B3=as.numeric(cybar$B3)
cybar$B4=as.numeric(cybar$B4)
cybar$B5=as.numeric(cybar$B5)
cybar$B6=as.numeric(cybar$B6)

rownames(cybar)<-Cybar_counts$X1

mat=vegdist(cybar)
#matrix mat is converted to dissimilarity matrix ('mat.dist') to be used with metaMDS
mat.cybar<-as.matrix(mat)         
# Renamed matrix
colnames(mat.cybar)<-rownames(mat.cybar)<-DATA3$Lago

#########################################################################
############################# Dalmation plot#############################
# Link to original code:
# http://www.ufz.de/index.php?en=38440
#########################################################################

# Reading files
dat<-read.delim("/home/bruno/Documentos/DOC/KOCH_PIPELINES/DAMATIONPLOT/Plots/Results_Greyscale.txt", row.names=2)
new.names<-gsub(".bmp","",row.names(dat))
row.names(dat)<-new.names

dat.ol<-dat[grep("overlap",row.names(dat)),]
dat.si<-dat[-c(grep("overlap",row.names(dat))),]

out.temp<-data.frame(t(combn(row.names(dat.si), 2)), t(combn(dat.si[,2], 2)))

out<-data.frame(out.temp, dat.ol[,2])
colnames(out)<-c("Scatter.1","Scatter.2","Pixel.Scatter.1","Pixel.Scatter.2","Pixel.overlap")
#jaccard index
#jacc.ind<- 2*out[,5]/(out[,3]+ out[,4])-1 # Original
jacc.ind<- out[,5]/((out[,3]+ out[,4])-out[,5]) # Modified
#jaccard distance
#jacc.dist<- (out[,5]-((out[,3]+out[,4])-out[,5]))/out[,5]  # Original
#jacc.dist<- 1- jacc.ind # Modified
jacc.dist<-  jacc.ind # Modified

#out.fin<-data.frame(out,Jaccard=jacc.dist) # Original
out.fin<-data.frame(out,Jaccard=jacc.dist) # Modified
size<-nrow(dat.si)
mat<-matrix(nrow=size, ncol=size)
j<-size-1
temp<-0
for (i in 1:j)  {
  mat[,i]<-c(rep(NA,size-length(out.fin[(temp+1):(temp+j),6])), out.fin[(temp+1):(temp+j),6])
  temp<-temp+j
  j<-j-1
}
diag(mat)<-0
colnames(mat)<-rownames(mat)<-rownames(dat.si)
mat.dalmation<-as.dist(mat) 


mat2=as.matrix(mat.dalmation)
mat.dp=mat2
mat.dp[upper.tri(mat.dp, diag = F)]=t(mat.dp)[upper.tri(mat.dp, diag = F)]


# Renamed matrix
colnames(mat.dp)<-rownames(mat.dp)<-DATA3$Lago

#########################################################################
############################# FlowFP ####################################
# Link to original code:
#http://www.ufz.de/index.php?en=38443
#########################################################################
library("flowFP")
mydata=samples$nodesample[[1]]

#Defines a model on sample 1 using the parameters FL1-H and SSC-H and 6 recursions.
model <- flowFPModel(mydata[[1]], parameters=c("FL1-H", "SSC-H"), nRecursions=6)

# Defining a flowSet from FCS (for FlowFP use)
pataFS=flowSet(mydata)
#Now, the model ist applied to all samples.
fp <- flowFP(pataFS, model)

#Displays the number of events per bin for all samples
# From here we can use vegan!
mat.fp0=counts(fp)
# Bray
mat.fp=vegdist(mat.fp0)
# Distance matrix
mat.fp=as.matrix(mat.fp)

# Renamed matrix
colnames(mat.fp)<-rownames(mat.fp)<-DATA3$Lago

#########################################################################
############################# FlowFDA####################################
#########################################################################
library("flowFDA")
#source("/home/bruno/Documentos/DOC/MAIN/Patagonia_OK/Props/MRM_parameters.R")
############################
################################################################################
### Function for sampling to equal nr. of cells in flowset
### Standard is minimum number of cells 
################################################################################

FCS.resample <- function(x, sample=0, replace=FALSE){
  library(easyGgplot2)
  sample_distr <- data.frame(counts=fsApply(x,FUN=function(x) nrow(x),use.exprs=TRUE))
  p1 <- ggplot2.histogram(data=sample_distr , xName='counts',
                          fill="white", color="black",
                          linetype="longdash",addMeanLine=TRUE, meanLineColor="red",
                          meanLineType="dashed", meanLineSize=1)+
    theme_bw() + labs(y="Frequency", title="Original count distribution")
  if(sample==0) sample <- min(fsApply(x=x,FUN=function(x) nrow(x),use.exprs=TRUE))
  ## Remove all .fcs files with less observations than the specified sample
  x <- x[fsApply(x=x,FUN=function(x) nrow(x),use.exprs=TRUE)>sample]
  for(i in 1:length(x)){
    exprs(x[[i]]) <- exprs(x[[i]])[sample(1:nrow(exprs(x[[i]])), sample, replace=replace),]
  }
  print(p1)
  cat(paste0("Your samples were randomly subsampled to ",sample," cells"))
  return (x)
}
### Beta-diversity (ranked based) using Non-metric Multidimensional Scaling (NMDS)
### x = flowBasis object from fingerprint (e.g., fingerprint)
### d = rounding factor for densities 
### n = number of replicates
### dist = choice of distance metric 
### k = number of MDS dimensions 
### iter = number of random starts in search of stable solution
beta.div.fcm <- function(x, d=3, n=1, dist="bray",k=2,iter=100,ord.type=c("NMDS","PCoA")){
  x <- x@basis/apply(x@basis, 1, max)
  require('vegan')
  input <- matrix(nrow=nrow(x)/n,ncol=ncol(x))
  j<-1
  if(n>1){
    for(i in seq(1,nrow(x),n)){
      if(n>1)
        input[j,] <- round(colMeans(x[(i:(i+(n-1))),]),d)
      j=j+1
    }
    rownames(input) <- rownames(x)[seq(1,nrow(x),n)]
    input.dist <- vegdist(input,method=dist)
    if(ord.type=="NMDS") mds.fbasis <- metaMDS(input.dist,autotransform=FALSE, k,trymax=iter)
    else mds.fbasis <- cmdscale(input.dist, k = 2, eig = TRUE, add = TRUE)
  }
  else{
    input.dist <- vegdist(x,method=dist)
    if(ord.type=="NMDS") mds.fbasis <- metaMDS(input.dist,autotransform=FALSE, k,trymax=iter)
    else mds.fbasis <- cmdscale(input.dist, k = 2, eig = TRUE, add = TRUE)
  }
  #return(mds.fbasis) # Original
  return(input.dist) # Modified
}

############################
set.seed(777)

flowData <- flowSet(samples[[2]][[1]])
sampleNames(flowData)<-samples$names
flowData_transformed <- flowData
param=c("FL1-H", "FL3-H","SSC-H","FSC-H")
flowData_transformed = flowData_transformed[,param]

### Randomly resample to the lowest sample size
flowData_transformed <- FCS.resample(flowData_transformed)
# V18 is removed in the ebove process

### Calculate fingerprint with bw = 0.01
fbasis <- flowBasis(flowData_transformed, param, nbin=128, 
                    bw=0.01,normalize=function(x) x)


### Beta-diversity assessment of fingerprint
beta.div <- beta.div.fcm(fbasis,n=1,ord.type="NMDS")

beta.div=as.matrix(beta.div)
dim(beta.div)
# Prepare to rename

# Renamed matrix
colnames(beta.div)<-rownames(beta.div)<-DATA3$Lago[-18]

#########################################################################
############################# PAIRWISE COMPARISONS#######################
#########################################################################
library("gplots")
# Mantel from all

out<-c(18,21,22,23,25) # Absent lakes in DGGE dataset and removed [18] througth FlowFDA process (resampling)
outprop<-c(20,21,22,24) # Absent lakes in DGGE in relation to to FlowFDA

# New matrices
DGGE.dm=as.matrix(vegdist(dgge[-18,-1]))
rownames(DGGE.dm)<-colnames(DGGE.dm)<-dgge$X1[-18]
DGGE.dm=as.dist(DGGE.dm)
CHIC.dm=as.dist(mat.chic2[-out,-out])
DALMATION.dm=as.dist(mat.dp[-out,-out])
CYBAR.dm=as.dist(mat.cybar[-out,-out])
FLOWFP.dm=as.dist(mat.fp[-out,-out])
FLOWDIV.dm=vegdist(m1[match(dgge$X1[-18], rownames(m1)),])
PROPS.dm=as.dist(beta.div[-outprop,-outprop])

# Heatmap function
myheatmap<-function(x){
  myfun <- function(x) hclust(as.dist(x), method = "ward.D")
  
  heatmap.2(x, hclustfun = myfun,  col=colorRampPalette(colors = c("white","azure","cyan2","blue4"))(100),
            dendrogram = "column", key=T, trace = "none",symm=T,
            lmat=rbind(c(0,3),c(2,1),c(0,4)), lwid=c(0.05,5), lhei = c(1,8,1.1), density.info = "none",
            key.title = NA,labCol=NA, key.xlab = NA, key.par=list(mar=c(2,0,0,8.3), cex=1),
            margins=c(.5,12.6),
            cexRow = 1.7,
            offsetRow = -.2)
}


myheatmap(as.matrix(CHIC.dm))
myheatmap(as.matrix(DALMATION.dm))
myheatmap(as.matrix(CYBAR.dm))
myheatmap(as.matrix(FLOWFP.dm))
myheatmap(as.matrix(FLOWDIV.dm))
myheatmap(as.matrix(PROPS.dm))
myheatmap(as.matrix(DGGE.dm))


# Mantel test

# Distance matrices
dms<-list(DGGE=DGGE.dm, CHIC=CHIC.dm, DALMATION=DALMATION.dm, CYBAR=CYBAR.dm, FLOWFP=FLOWFP.dm, PROPS=PROPS.dm, FLOWDIV=FLOWDIV.dm)
lapply(dms, function(x)dim(as.matrix(x)))

mantel.t=list()
length(mantel.t)=7
names(mantel.t)<-c("DGGE", "CHIC", "DALMATION", "CYBAR", "FLOWFP", "PROPS", "FLOWDIV")

for (i in 1:7)mantel.t[[i]]=lapply(dms, function(x)(mantel(x, dms[[i]])))

mantel.t

mantel.t2=lapply(mantel.t, function(x)lapply(x, function(y)y$signif))
mantel.t3=lapply(mantel.t, function(x)lapply(x, function(y)y$statistic))


output1 <- matrix(unlist(mantel.t2), ncol = 7, byrow = TRUE)
colnames(output1)<-rownames(output1)<-names(mantel.t)
print("Significance")
output1

output2 <- matrix(unlist(mantel.t3), ncol = 7, byrow = TRUE)
colnames(output2)<-rownames(output2)<-names(mantel.t)
print("Statistics")
output2

