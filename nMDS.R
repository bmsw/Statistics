
#Codes for nMDS

colorcode=substring(id.flowdiv_clean, 1, 2)

colorcode[colorcode%in%c("01", "05", "16", "31")]= "cadetblue1"

colorcode[colorcode%in%c("21", "25", "26", "33")]= "chocolate1"

colorcode[colorcode%in%c("06", "08", "14", "20")]= "darkorchid4"

colorcode[colorcode%in%c("18", "19")]= "deepskyblue"

colorcode[colorcode%in%c("15", "17")]= "magenta"

colorcode[colorcode%in%c("23", "24")]= "green4"

colorcode[colorcode%in%c("22", "09")]= "hotpink4"

colorcode[colorcode%in%c("28", "30")]= "yellow2"

colorcode[colorcode%in%c("12", "11", "13", "29", "03", "32", "02", "27", "07")]= "black"


#nMDS
meta1=metaMDS(m1, autotransform = T)

meta2=metaMDS(m1, autotransform = F)


plot(meta1, type="n", main="flowDiv nMDS")
text(meta1$points,substring(id.flowdiv_clean, 1, 2), col=colorcode, pch=19)
plot(meta2, type="n", main="flowDiv nMDS")
text(meta2$points,substring(id.flowdiv_clean, 1, 2), col=colorcode, pch=19)
