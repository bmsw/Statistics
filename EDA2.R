library(lattice)
library(directlabels)
library(ggplot2)

#Seleção de variável
data0=dados_
var1="OCS.n"
#var1="FAL.n"
var2="ano.1"

# Boxplot
ggplot(aes_string(y=var1, x='1'), data=data0)+
  geom_boxplot()+
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  labs(y=toupper(var1), x="")


# Distribution
descdist(data0[,var1], discrete = TRUE)
# Mosaic plot
tb2=table(data0[,var2], factor(data0[,var1], levels=0:max(data0[,var1])))
plot(tb2, main="Mosaic Plot", col=alpha("blue", 0.5))
# Barplot all
tb1=table(factor(data0[,var1], levels=0:max(data0[,var1])))
d1<-data.frame(tb1)
d1$Var1<-as.numeric(as.character(d1$Var1))
g1<-ggplot(aes(x=Var1, y=Freq), data=d1,stat ="identity", position="dodge")+
  geom_bar(stat ="identity", position="dodge", col="blue", fill="green",alpha = .5)+
  geom_line(col="red", size=.5)+
  geom_point()+
  scale_x_continuous(breaks = c(0:(length(d1$Var1)-1)))+
  labs(y="Contagem", x="Números de eventos")
g1
  
# Barplot by factor
tb3<-data.frame(tb2)
#Prop
a=prop.table(tb2)
b=data.frame(a)
# Grid
# Absoluto
xyplot(Freq~Var2|Var1, data=tb3, type=c("h", "p"), pch=16, lwd=2, cex=.75, xlab="Números de eventos", ylab="Contagem")
#Proporção
xyplot(Freq~Var2|Var1, data=b, type=c("h", "p"), pch=16, lwd=2, cex=.75, xlab="Números de eventos", ylab="Frequência")
# All together
#Absoluto
plot1<-xyplot(Freq~Var2, groups=Var1, data=tb3, type="b", pch=16, lwd=2, cex=.75, xlab="Números de eventos", ylab="Contagem")
direct.label(plot1, list("top.points", cex=1, dl.trans(y=y+.1)))
#Proporção
plot1<-xyplot(Freq~Var2, groups=Var1, data=b, type="b", pch=16, lwd=2, cex=.75, xlab="Números de eventos", ylab="Frequência")
direct.label(plot1, list("top.points", cex=1, dl.trans(y=y+.1)))


