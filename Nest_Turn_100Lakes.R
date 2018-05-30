data_bray<-data.frame(Nestedness=nest[1,]/bray[1,], Turnover=turn[1,]/bray[1,], Distancia=distancias[1,], Lon=dataall$Lon, Bray=bray[1,], Lat=dataall$Lat)
data_bray2<-data_bray[-1,]

# 
# dataall.cart<-geoXY(dataall$Lat, dataall$Lon, unit=1000)
# distancias<-as.matrix(dist(dataall.cart, "euclidian"))
# d0<-distancias[upper.tri(distancias)]

toy<-data_bray2
toy2<-stack(toy[,c(1:2)])
toy3<-data.frame(toy2, dis=rep(toy$Distancia,2))
ggplot(toy3, aes(x=dis, y=values, fill=ind)) + geom_area()

ggplot(toy3, aes(x=dis, y=values, fill=ind)) + geom_line()


n1<-nest/bray
t1<-turn/bray

dgeo<-distancias[upper.tri(distancias)]
dnest<-n1[upper.tri(n1)]
dturn<-t1[upper.tri(t1)]


d1<-data.frame(dgeo, dnest, dturn)
d2<-stack(d1[,c(2:3)])
d3<-data.frame(d2, dis=rep(d1$dgeo,2))

ggplot(d3, aes(x=dis, y=values, fill=ind))  + geom_area()+
  geom_smooth(aes(group=ind, linetype=ind), size=.8, color="gray30", se=F, method = "loess")+
  scale_x_continuous(limits = c(0, 150))
