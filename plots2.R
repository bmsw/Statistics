#library(viridis)
#library(cowplot)


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
    axis.title=element_text(size=30),
    # axis.line = element_line(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  theme(legend.position="none")
###################################
g2<- ggplot(tubed[[i]])+
  aes_string(x="SSC.H", y="FL1.H")+
  geom_hex(bins = 100, na.rm = T) +
  scale_fill_gradientn("", colours = rev(gray.colors(3, end = 4/6)))+
  scale_y_continuous(breaks=breaks[,"FL1.H"],  limits= c(myseqs$max.min.table["FL1.H", "minimum"], myseqs$max.min.table["FL1.H", "maximum"])) +
  scale_x_continuous(breaks=breaks[,"SSC.H"],  limits= c(myseqs$max.min.table["SSC.H", "minimum"], myseqs$max.min.table["SSC.H", "maximum"])) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        axis.line = element_blank(),
        axis.title=element_text(size=30),
        # axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  theme(legend.position="none")+
  geom_point(data=clusters2, aes_string(y="FL1.H", x="SSC.H"), colour=clusters2$clust, alpha=.5, shape=15, size=1.83)


g3<-ggdraw() +
  draw_plot(g2, 0.4, 0, 0.6, 0.6) +
  draw_plot(g1, 0, 0.5, .5, .5) 

ggsave(paste(i,".png" ), g3, path=getwd(), device="png", units="in", width = 9.8, height = 9.8, dpi=300)
}

#ggsave("teste.png",g3, path=getwd(), device="png", units="in", width = 9.8, height = 9.8, dpi=300)