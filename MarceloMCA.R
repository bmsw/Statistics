library("FactoMineR")
library("factoextra")

marcelo_data=read.xlsx("/home/brunomsw/Documentos/DOUTORADO/Projetos[Parcerias]/Marcelo/tab_especie_mes_arte_bruta.xlsx",2)

marcelo_data=marcelo_data[,-c(2,5)]
#marcelo_data=marcelo_data[,c(2,5)]
summary(marcelo_data)

  for (i in 1:5) {
    library(forcats)
    
    df1=data.frame(table(marcelo_data[,i]))
    
    p <- ggplot(df1, aes(x = Var1, y = Freq)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = Freq), vjust = 1.5, colour = "white")+
      labs(title=colnames(marcelo_data)[i], y="Count")
    
    
    print(p + aes(x = fct_reorder(Var1, Freq))+xlab)
    
  }


res.mca <- MCA(marcelo_data, graph = FALSE)

print(res.mca)

eig.val <- get_eigenvalue(res.mca)

fviz_screeplot(res.mca, addlabels = TRUE, ylim = c(0, 15), title="")

#fviz_mca_biplot(res.mca, 
                #repel = TRUE, # Avoid text overlapping (slow if many point)
                #ggtheme = theme_minimal())
fviz_mca_var(res.mca, 
repel = TRUE, # Avoid text overlapping (slow if many point)
ggtheme = theme_minimal(), col.var = "black", title="")

fviz_cos2(res.mca, choice = "var", axes = 1:2)

#plot(res.mca, invisible = "ind", col.var = "black")

# Cos2 of variable categories on Dim.1 and Dim.2
fviz_cos2(res.mca, choice = "var", axes = 1:2)

# Contributions of rows to dimension 1
fviz_contrib(res.mca, choice = "var", axes = 1, top = 15)
# Contributions of rows to dimension 2
fviz_contrib(res.mca, choice = "var", axes = 2, top = 15)

