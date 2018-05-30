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


 
  myworkspaces="/home/bruno/Documentos/DOC/MAIN/Patagonia_OK/METADATA/flowDiv patagonia/Romina_Bruno_1_10.wsp"
  gate.name="Bact tot"; ialpha="invsimpson"; ibeta="bray"; do.plot=FALSE; static=FALSE; dilutions=DATA$`x (flowDiv)`;  pmax=0.05; env=NULL; use.beads=FALSE; beads=NULL; transform.hell=FALSE; dimension=20; psize=3; autotrans=FALSE
    
  
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

  
  expressions=lapply(unlist(nnn$nodesample), function(x)exprs(x)[,selection])
  tubed=lapply(expressions, function(x)data.frame(x))
  
  breaks=apply(myseqs$cuts, 2, function(x){
    bind=cbind(x[1:length(x)-1], x[2:length(x)]) # Bindinfg for processing on next statge
    brk=apply(bind, 1, mean)
    return(brk)
  }
  )

  tubed=tubed
  names=nnn$names
  NBINS=nbins
  selection=selection
  myseqs=myseqs
  breaks=breaks
  psize=psize
  dimension=dimension 
  
  breaks=data.frame(breaks)
  myseqs=lapply(myseqs, function(x)data.frame(x))
  rownames(myseqs$max.min.table)<-gsub( "-", ".", rownames(myseqs$max.min.table))
 
  
  ggplot(tubed[[6]])+
    aes(x=SSC.H, y=FL1.H)+
    geom_hex(bins = 100, na.rm = T) +
    scale_fill_gradientn("", colours = rev(gray.colors(3, end = 4/6)))+
    scale_y_continuous(breaks=breaks[,"FL1.H"],  limits= c(myseqs$max.min.table["FL1.H", "minimum"], myseqs$max.min.table["FL1.H", "maximum"])) +
    scale_x_continuous(breaks=breaks[,"SSC.H"],  limits= c(myseqs$max.min.table["SSC.H", "minimum"], myseqs$max.min.table["SSC.H", "maximum"])) +
    theme_bw() +
    theme(#panel.border = element_blank(), 
          axis.line = element_blank(),
          axis.title=element_text(size=40),
          # axis.line = element_line(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    theme(legend.position="none")
  
  
################################
  

    
################################  
######SAving all
  
  
  for (i in 1:31){
    
    g1<-ggplot(tubed[[i]])+
      aes(x=SSC.H, y=FL1.H)+
      geom_hex(bins = 100, na.rm = T) +
      scale_fill_gradientn("", colours = rev(gray.colors(3, end = 4/6)))+
      scale_y_continuous(breaks=breaks[,"FL1.H"],  limits= c(myseqs$max.min.table["FL1.H", "minimum"], myseqs$max.min.table["FL1.H", "maximum"])) +
      scale_x_continuous(breaks=breaks[,"SSC.H"],  limits= c(myseqs$max.min.table["SSC.H", "minimum"], myseqs$max.min.table["SSC.H", "maximum"])) +
      theme_bw() +
      theme(#panel.border = element_blank(), 
        axis.line = element_blank(),
        axis.title=element_text(size=40),
        # axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
      theme(legend.position="none")
    # 100 = 1 polegada no Rstudio
    ggsave(paste("r",i,".png" ), g1, path=getwd(), device="png", units="in", width = 5.96, height = 5.9, dpi=300)
  }
  