library(Amelia)
library(gvlma)
library(beanplot)
library(ggplot2)
library(sm)
# Importar os dados
setwd("/home/bruno/Documentos/DOC/Parcerias/Marcelo/artigoeplanilhadedadosr_porosus/Art03")
mar_data=read.csv("plan_data_brutoCSV.csv", sep="\t", header=T, dec=",")
missmap(mar_data, legend = F, main = "")

# Limpeza (NA's e Outliers)
#mar_data=mar_data[-c(51,107),complete.cases(t(mar_data))]
mar_data=mar_data[-c(51,107),] # Outliers
#NA's para Pevis
mar_data=mar_data[!is.na(mar_data$Pevis),]
missmap(mar_data, legend = F, main = "")


attach(mar_data)

mar_data$logCT<-log10(CT)
mar_datalogCP<-log10(CP)

logCT=mar_data$logCT
logCP=mar_datalogCP
#Visual check

plot(logCT, logCP, type="n")
text(logCT, logCP, labels = as.character(c(1:length(CP))))

plot(log10(Pevis), logCP, type="n")
text(log10(Pevis), logCP, labels = as.character(c(1:length(CP))))

# Ancova Não Paramétrica
# 
sm.ancova(logCP,log10(Pevis), Sexo, model="equal") # Significativo!
sm.ancova(logCP,log10(Pevis), Observacao, model="equal") # 

sm.ancova(log10(Pevis),log10(Pevis), Observacao, model="equal") # 


summary(lm(log10(Pevis)~logCP+Sexo))
sm.ancova(CP, Pevis,  Sexo, model="equal", h=3) # 



sm.ancova(CP, Pevis, Observacao, model="equal", h=3) # 




# msexo<-lm(log10(Pevis)~logCP+Sexo)
# gvlma(msexo)
# summary(msexo)
# 
# mobs<-lm(log10(Pevis)~logCP+Observacao)
# gvlma(mobs)
# summary(mobs)
# 
# #Beter way
# msexo2<-aov(logCT~logCP*Sexo)
# msexo3<-aov(logCT~logCP+Sexo)
# 
# summary(msexo2)
# summary(msexo3)
# 
# 
# anova(msexo2, msexo3)
# 
# mobs2<-aov(logCT~logCP*Observacao)
# mobs3<-aov(logCT~logCP+Observacao)
# 
# summary(mobs2)
# summary(mobs3)
# 
# anova(mobs2, mobs3)

# Ggplot

ggplot(mar_data, aes(x=CT, y=CP, shape=Observacao, colour=Observacao))+geom_point()+geom_smooth(method="lm", fill=NA)
ggplot(mar_data, aes(x=CT, y=CP, shape=Sexo, colour=Sexo))+geom_point()+geom_smooth(method="lm", fill=NA)

#Chi quadrado
table01<-table(Observacao, Sexo)
plot(table01, main="", cex=1.2)
chisq.test(table01)

# Kruskall-Wallis
kruskal.test(Pevis, Sexo)
kruskal.test(Pevis, Observacao)
dunn.test(Pevis, Observacao, method = "Bonferroni")


kruskal.test(CT, Sexo)
kruskal.test(CT, Observacao)

dunn.test(CT, Observacao, method = "Bonferroni")

# Beanplots

beanplot(CT~Sexo, what = c(1,1,1,0), log = "")
beanplot(CT~Observacao, what = c(1,1,1,0), log = "")
