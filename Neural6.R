#Fósforo total variou de 0.01 a 2.63 μgl-1
# clorofila a de 0.15 μgl-1 a 504.9 μgl-1
mydata2<-subset(dataall, select = c(colunas, 49, 33, 20, 34))
with(mydata2, tapply(ph, cat_pH, summary))
with(mydata2, tapply(Ptot, Cat_PT, summary))
with(mydata2, tapply(LOGCLA, Cat.cla, summary))
with(mydata2, tapply(LOGCLA, Cat_cla, summary))
with(mydata2, tapply(cla, Cat.cla, summary))
with(mydata2, tapply(cla, Cat_cla, summary))

# for chla
x=dataall$cla
x=replace(x, x<2.6, 1)
x=replace(x, x>=2.6&x<20, 2)
x=replace(x, x>=20&x<56, 3)
x=replace(x, x>=56, 4)
dataall$cla.carlson<-factor(x)
# for P
x=dataall$Ptot
x=replace(x, x<12, 1)
x=replace(x, x>=12&x<24, 2)
x=replace(x, x>=24&x<96, 3)
x=replace(x, x>=96, 4)
dataall$ptot.carlson<-factor(x)
# for secchi
x=dataall$Secc..m.
x=replace(x, x<.5, 4)
x=replace(x, x>=.5&x<2, 3)
x=replace(x, x>=2&x<4, 2)
x=replace(x, x>=4, 1)
dataall$secchi.carlson<-factor(x)
