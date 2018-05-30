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
DATA$`TROPHIC STATE`=revalue(DATA3$`TROPHIC STATE`, c("EUT"="MESO_EUT", "MESO"="MESO_EUT", "OLIGO"="OLIGO"))
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
# 3 dimensions, 67 bins per dimesnion
#pata.fd=flowDiv("/home/bruno/Documentos/DOC/MAIN/Patagonia_OK/METADATA/flowDiv patagonia/Romina_Bruno_1_10.wsp", "Bact tot", do.plot = F, static = F, dilutions=DATA$`x (flowDiv)`, transform = T, autotrans=T,psize = 1.08)


myworkspaces="/home/bruno/Documentos/DOC/MAIN/Patagonia_OK/METADATA/flowDiv patagonia/Romina_Bruno_1_10.wsp"
  gate.name="Bact tot"
  do.plot = T
  static = F
  dilutions=DATA$`x (flowDiv)`
  transform = T
  autotrans=T
  psize = 1.08
  ialpha="invsimpson"; ibeta="bray"; pmax=0.05; use.beads=FALSE; beads=NULL; transform.hell=FALSE; dimension=20
  
  classes=sapply(myworkspaces, function(x) return(class(x)), simplify = T)
  gsets=which(classes=="GatingSet")
  if(length(gsets)>0){
    wksp1=myworkspaces[-gsets]
    gtsets=myworkspaces[gsets]
    wksp2<-opc(wksp1)
    wksp<-append(wksp2, gtsets)
    
  } else (wksp=opc(myworkspaces))
  nnn=nn(wksp, nod=gate.name,use.beads=use.beads, nod2=beads)
  
  ops<-unique(unlist(lapply(unlist(nnn$nodesample), colnames)))
  selection<-select.list(ops, c(1:length(ops)), multiple = T, title="Please select channels for use")
  
  fixed=NULL
  if (static==TRUE){
    fixeds=lapply(selection, function(x){
      message(paste("Please enter minimum value for", x))
      min=scan(nmax=1, quiet=T)
      
      message(paste("Please enter maximum value for", x))
      max=scan(nmax=1, quiet=T)
      
      return(c(minimum=min, maximum=max))
    } )
    
    fixed=do.call(rbind, fixeds)
    rownames(fixed)<-selection
  }
  
  
  binss=lapply(nnn$nodesample, function(x) lapply(x, function(y)lapply(selection, function(x)Freedman.Diaconis(exprs(y)[,x]))))
  suggested.bins=round(median(unlist(binss)))
  
  message(paste("Suggested number of bins:", suggested.bins, "\n"), "How many bins do you want to use?")
  nbins=scan(nmax=1, quiet=T)
  
  myseqs <- seq_fun(nnn$nodesample, selection, nbins, fixed)
  
  hits<-lapply(nnn$nodesample, function(x) lapply(x, function(y)sapply(selection, function(x){
    cuts<-cut(exprs(y)[,x], myseqs$cuts[,x],include.lowest = T, labels=F)
    return(cuts)
  }, simplify=F)))
  hits=unlist(hits, recursive=FALSE)
  bin.levels=lapply(hits, function(x)lapply(x, function(y)factor(as.numeric(y), levels=1:nbins)))
  tables=lapply(bin.levels, function(x)table(x))
  
  if (!is.null(dilutions))tables=mapply(function(x,y)x*y, tables, dilutions, SIMPLIFY = F)
  
  matrices<-do.call(mapply, c(cbind, tables))
  
  colnames(matrices)<-as.character(c(1:nbins^length(selection)))
  rownames(matrices)<-unlist(nnn$names, recursive = F)
  
  
  if (do.plot==TRUE){
    
    
    expressions=lapply(unlist(nnn$nodesample), function(x)exprs(x)[,selection])
    tubed=lapply(expressions, function(x)data.frame(x))
    
    breaks=apply(myseqs$cuts, 2, function(x){
      bind=cbind(x[1:length(x)-1], x[2:length(x)]) # Bindinfg for processing on next statge
      brk=apply(bind, 1, mean)
      return(brk)
    }
    )
    
    set.seed(1)
    nmds=metaMDS(matrices,autotransform = autotrans, try=100)
    specs=na.omit(nmds$species)
    
    
    ccKM=cascadeKM(specs, 2 ,10)
    ngroups=names(which(ccKM$results["SSE",]==min(ccKM$results["SSE",])))
    ngroups=gsub(" groups","", ngroups)
    message(paste("Suggested number of clusters:", ngroups, "\n"), "Please enter the number of clusters to use:")
    nclusters=scan(nmax=1, quiet=T)
    
    km=kmeans(specs, nclusters)
    
    cbins=as.numeric(names(km$cluster))
    bb=sapply(cbins, function(x)binpos(tables[[1]], x), simplify = F)
    bb=do.call(rbind, bb)
    clusters=lapply(selection, function(x)breaks[bb[,x],x])
    clusters=t(do.call(mapply, c(rbind, clusters)))
    colnames(clusters)<-selection
    clusters<-data.frame(clusters)
    
    
    invisible(capture.output(mapply(graphs, tubed=tubed, names=nnn$names, MoreArgs = list(NBINS=nbins,clusterscols=km$cluster, selection=selection, myseqs=myseqs, breaks=breaks, clusters=clusters, psize=psize, dimension=dimension))))
    
    plot(nmds, type="n")
    points(specs, pch=19, col=adjustcolor(km$cluster, alpha.f = 0.3))
    text(nmds$points, labels = c(1:nrow(matrices)))
    if (!is.null(env)){
      fit <- envfit(nmds, env, perm = 999)
      plot(fit, p.max = pmax, col = "red")
    }
    
    
    
    
    
  }
  
  alpha<-apply(matrices, 1, function(x)diversity(x, index=ialpha))
  
  pielou<-apply(matrices, 1, function(x)diversity(x)/log(vegan::specnumber(x)))
  
  if(transform.hell){
    beta<-vegdist(decostand(matrices, "hell"), method=ibeta)
  }
  
  else {
    if(ibeta=="bray"){
      beta<-bray.part(matrices)
    }
    else {
      beta<-as.matrix(vegdist(matrices, method=ibeta))
    }
  }
  

