library(gvlma)
library(ggplot2)
load("~/Data100Lakes.RData")
# Usando o log do pCO2
dataall$Logpco2<-log10(dataall$pco2..µatm.)
#dataall$rich<-log10(dataall$rich) # Mudar na linha 121 para "Log10 Richness"
# Nova variavel
#names=c("a","d","d","e","f",colnames(mydata))
######### ALFA
mydata=dataall # New dataset (to not override!)
mydata=mydata[,c(9,13:37,48:52)]

###################################
library(fitdistrplus)
  with(mydata, descdist(pielou))
  with(mydata, descdist(rich))
  with(mydata, descdist(shannon))

with(mydata, plot(pielou, shannon))
with(mydata, plot(rich, shannon))
with(mydata, plot(rich, pielou))

####################################
# Using justo 10:20 >>> Change it!

# for (i in 1:30){
#   bc=boxCox(lm(mydata[,i]~mydata$pielou), family="yjPower")
#   best=with(bc, x[which.max(y)])
#   mydata[,ncol(mydata)+1]<-yjPower(mydata[,i], best)
#   colnames(mydata)[ncol(mydata)]<-paste("TransP", colnames(mydata)[i], sep="_")
# }

##################

## Função para verificar a adequação dos modelos lineares (gvlma OK pra todos os testes)
gvbruno<-function(a){
  soma=sum(a$GlobalTest$GlobalStat4$Decision, a$GlobalTest$DirectionalStat1$Decision, a$GlobalTest$DirectionalStat2$Decision, a$GlobalTest$DirectionalStat3$Decision, a$GlobalTest$DirectionalStat4$Decision)
  return(soma==0)}

# Verificando gvlma
gr=c(); gs=c(); gp=c(); gsi=c(); ginv=c()
r=c(); p=c(); sh=c(); si=c(); inv=c()
for(i in 1:31) {
  # Verificando gvlma
  gr[i]=gvbruno(gvlma(lm(rich~mydata[,i])))
  gp[i]=gvbruno(gvlma(lm(pielou~mydata[,i])))
  gs[i]=gvbruno(gvlma(lm(shannon~mydata[,i])))
  gsi[i]=gvbruno(gvlma(lm(simpson~mydata[,i])))
  ginv[i]=gvbruno(gvlma(lm(invsimpson~mydata[,i])))
  # Verificando Correlações
  r[i]=with(dataall, cor.test(rich, mydata[,i]))$p.value
  p[i]=with(dataall, cor.test(pielou, mydata[,i]))$p.value
  sh[i]=with(dataall, cor.test(shannon, mydata[,i]))$p.value
  si[i]=with(dataall, cor.test(simpson, mydata[,i]))$p.value
  inv[i]=with(dataall, cor.test(invsimpson, mydata[,i]))$p.value
  
  
}

#gvlma(lm(dataall$shannon~dataall$PRECIPITAÇÃO.ANUAL)) # ???
# Verificando Correlações


# Dados das correlações <= 0.05 e gvlma's correspondentes
correls=cbind(cbind(p, r, sh, si, inv)<=0.05, gp, gr, gs,gsi, ginv, colnames(mydata))
correls<-data.frame(correls)
# Organizando tipos de dados para booleanos
for (i in 1:(ncol(correls)-1))correls[,i]<-as.logical(correls[,i])

# Snooping: Verificando modelos significativos && válidos (segundo gvlma)
snoop1=cbind(p=(correls$p)&(correls$gp), r=(correls$r)&(correls$gr), sh=(correls$sh)&(correls$gs), si=(correls$si)&(correls$gsi), inv=(correls$inv)&(correls$ginv))
snoop2=data.frame(snoop1, colnames(mydata))
##
snoop3=snoop2[!apply(snoop1,1, sum)==0,]

# Pielou
snoop4<-snoop3[snoop3[,1]==TRUE,]
snoop4[,6]<-as.character(snoop4[,6])
snoop4$labels<-c("CO2 Pressure", "Secchi", "430nm Absorbance", "Log10 COD", "Dissolved P", "Chla", "Atomic P", "Log10 Chla", "Log10 CO2", "Log10 Chla", "Log10 pCO2")
snoop4

for (i in 1:nrow(snoop4)){
  print(summary(lm(dataall$pielou~dataall[,snoop4[i,6]])))
  #############
  m<-lm(dataall$pielou~dataall[,snoop4[i,6]])
  eq <- substitute(italic(y) == a*sinal*b %.% italic(x)*~~~~~~"("~italic(r)^2~"="~r2~")", 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        sinal = ifelse(coef(m)[2]<0, "", "+"),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  regeq<-as.character(as.expression(eq))
  ####################### 
  g1<-ggplot(aes_string(snoop4[i,6], "pielou"), data=dataall)+geom_point()+geom_smooth(method='lm')+
    ylab("Pielou Index")+xlab(snoop4[i,7])+
    ggtitle(snoop4[i,7], subtitle =  as.formula(regeq))+
    theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))
  
  
  plot(g1)
  
}


# Richness
snoop4<-snoop3[snoop3[,2]==TRUE,]
snoop4[,6]<-as.character(snoop4[,6])
snoop4$labels<-c("Log10 Chla","Log10 Chla")

for (i in 1:nrow(snoop4)){
  print(summary(lm(dataall$rich~dataall[,snoop4[i,6]])))
  #############
  m<-lm(dataall$rich~dataall[,snoop4[i,6]])
  eq <- substitute(italic(y) == a*sinal*b %.% italic(x)*~~~~~~"("~italic(r)^2~"="~r2~")", 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        sinal = ifelse(coef(m)[2]<0, "", "+"),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  regeq<-as.character(as.expression(eq))
  ####################### 
  g1<-ggplot(aes_string(snoop4[i,6], "rich"), data=dataall)+geom_point()+geom_smooth(method='lm')+
    ylab("Richness")+xlab(snoop4[i,7])+
    ggtitle(snoop4[i,7], subtitle =  as.formula(regeq))+
    theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))
  
  
  plot(g1)
  
}


