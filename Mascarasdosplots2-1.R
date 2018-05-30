

tubed=tubed
names=nnn$names
NBINS=nbins
clusterscols=brewer.pal(8, "Dark2")[km$cluster]
selection=selection
myseqs=myseqs
breaks=breaks
clusters=clusters
psize=psize
dimension=dimension

#Settings for saving plot with those configs:
# Heigth=504 Weigth=494


  
  # newframe=expand.grid(selection, selection)
  # newframe2=as.matrix(newframe)
  # newframe2= gsub("-",".", newframe2)
  # 
  breaks=data.frame(breaks)
  myseqs=lapply(myseqs, function(x)data.frame(x))
  rownames(myseqs$max.min.table)<-gsub( "-", ".", rownames(myseqs$max.min.table))
  clusters2=cbind(clusters, clust=clusterscols)
  
  #clusters2=clusters2[!duplicated(clusters2[,c("SSC.H","FL1.H")]),]
   
    
  # 11 
  #28
  #6
  #30

  
  # Ficam: 11 (P13), 6 (P7)
      
    g1<- ggplot(tubed[[6]])+
        aes_string(x="SSC.H", y="FL1.H")+
        geom_hex(bins = 100, na.rm = T) +
        scale_fill_gradientn("", colours = rev(gray.colors(3, end = 4/6)))+
        scale_y_continuous(breaks=breaks[,"FL1.H"],  limits= c(myseqs$max.min.table["FL1.H", "minimum"], myseqs$max.min.table["FL1.H", "maximum"])) +
        scale_x_continuous(breaks=breaks[,"SSC.H"],  limits= c(myseqs$max.min.table["SSC.H", "minimum"], myseqs$max.min.table["SSC.H", "maximum"])) +
        theme_bw() +
        theme(panel.border = element_blank(), 
              axis.line = element_blank(),
              axis.title=element_text(size=40),
             # axis.line = element_line(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())+
        theme(legend.position="none")+
        geom_point(data=clusters2, aes_string(y="FL1.H", x="SSC.H"), colour=clusters2$clust, alpha=.5, shape=15, size=1.83)
     
     
    # Só  maks!
     
     
     ggplot(tubed[[6]])+
scale_fill_gradientn("", colours = rev(gray.colors(3, end = 4/6)))+
       scale_y_continuous(breaks=breaks[,"FL1.H"],  limits= c(myseqs$max.min.table["FL1.H", "minimum"], myseqs$max.min.table["FL1.H", "maximum"])) +
       scale_x_continuous(breaks=breaks[,"SSC.H"],  limits= c(myseqs$max.min.table["SSC.H", "minimum"], myseqs$max.min.table["SSC.H", "maximum"])) +
       theme_bw() +
       theme(panel.border = element_blank(), 
             axis.title=element_text(size=40),
             axis.line = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(), 
             axis.text.x=element_blank(),
             axis.ticks.x=element_blank(),
             axis.text.y=element_blank(),
             axis.ticks.y=element_blank())+
       theme(legend.position="none")+
       labs(x = "SSC-H", y = "FL1-H")+
       geom_point(data=clusters2, aes_string(y="FL1.H", x="SSC.H"), colour=clusters2$clust, alpha=.5, shape=15, size=1.83)
     
    
     plot(nmds, type="n")
     points(specs, pch=15, cex=1.5, col=adjustcolor(clusterscols, alpha.f = 0.2))
     text(nmds$points, cex=1.5, labels = c(1:nrow(matrices)))
     
    
     
     
     
     
     
     
     
     
     
     ##### para salvar todos os graficos
     
     for (i in 1:31){
       
       g1<-ggplot(tubed[[i]])+
         aes_string(x="SSC.H", y="FL1.H")+
         geom_hex(bins = 100, na.rm = T) +
         scale_fill_gradientn("", colours = rev(gray.colors(3, end = 4/6)))+
         scale_y_continuous(breaks=breaks[,"FL1.H"],  limits= c(myseqs$max.min.table["FL1.H", "minimum"], myseqs$max.min.table["FL1.H", "maximum"])) +
         scale_x_continuous(breaks=breaks[,"SSC.H"],  limits= c(myseqs$max.min.table["SSC.H", "minimum"], myseqs$max.min.table["SSC.H", "maximum"])) +
         theme_bw() +
         theme(panel.border = element_blank(), 
               axis.line = element_blank(),
               axis.title=element_text(size=40),
               # axis.line = element_line(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), 
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               axis.text.y=element_blank(),
               axis.ticks.y=element_blank())+
         theme(legend.position="none")+
         geom_point(data=clusters2, aes_string(y="FL1.H", x="SSC.H"), colour=clusters2$clust, alpha=.5, shape=15, size=1.83)
       
       # 100 = 1 polegada no Rstudio
       ggsave(paste("L.",i,"_", DATA$Lago[i], ".png" ), g1, path=getwd(), device="png", units="in", width = 5.96, height = 5.9, dpi=300)
     }
     