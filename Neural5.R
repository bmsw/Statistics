# Usado com o script neural 3 e 4 para criação das matrizes de confusão
library(caret)

colunas= c(38, 39, 42, 43, 45)
preditores=c(1,3)
# 1 e 3
(mydata<-dataall[,c(preditores, colunas)])


neuralbruno<-function(mydata, layers=5){
  mydata.<-mydata
  #mydata.<-teste
  # turnin into factor
  mydata.[,(length(preditores)+1)]<-factor(mydata.[,(length(preditores)+1)])
  # training set
  train<-cbind(mydata.[,-(length(preditores)+1)], class.ind(mydata.[,(length(preditores)+1)]))
  #renamiing
  nlabels<-ncol(train)-length(preditores)
  
  if(nlabels==3){
    colnames(train) <- c(names(mydata.)[1:length(preditores)],"l1","l2","l3")
  }
  if(nlabels==4){
    colnames(train) <- c(names(mydata.)[1:length(preditores)],"l1","l2","l3", "l4")
  }
  if(nlabels==5){
    colnames(train) <- c(names(mydata.)[1:length(preditores)],"l1","l2","l3", "l4", "l5")
  }
  
  #scaling
  scl <- function(x){ (x - min(x))/(max(x) - min(x)) }
  if(length(preditores)==1){
    train[, 1] <- scl(train[, 1])
  } else{
    train[, c(1:length(preditores))] <- data.frame(lapply(train[, c(1:length(preditores))], scl))
  }
  
  # renaming and formula to neuralnet
  n <- colnames(train)
  
  if(nlabels==3){
    f <- as.formula(paste("l1 + l2 + l3 ~", paste(n[!n %in% c("l1","l2","l3")], collapse = " + ")))
  }
  
  if(nlabels==4){
    f <- as.formula(paste("l1 + l2 + l3 + l4 ~", paste(n[!n %in% c("l1","l2","l3", "l4")], collapse = " + ")))
  }
  
  if(nlabels==5){
    f <- as.formula(paste("l1 + l2 + l3 + l4 + l5 ~", paste(n[!n %in% c("l1","l2","l3", "l4", "l5")], collapse = " + ")))
  }
  
  nn <- neuralnet(f,
                  data = train,
                  hidden = layers,
                  act.fct = "logistic",
                  linear.output = FALSE,
                  lifesign = "minimal", threshold = .05)
  
  #plot(nn)
  # Compute predictions
  pr.nn <- compute(nn, train[, c(1:length(preditores))])
  
  # Extract results
  pr.nn_ <- pr.nn$net.result
  #head(pr.nn_)
  
  #original_values <- max.col(train[, nlabels:ncol(train)])
  original_values <- max.col(train[, (length(preditores)+1):ncol(train)])
  pr.nn_2 <- max.col(pr.nn_)
  media<-mean(pr.nn_2 == original_values)
  
  ######################################################
 
  return(data.frame(original=factor(original_values), predicted=factor(pr.nn_2)))
}


#"Cat_cla"
teste<-mydata[,c(1,2,3)]
table(teste[,3])
cm<-neuralbruno(teste, layers = 4)
confusionMatrix(cm$predicted, cm$original)
#"Cat_PT"
teste<-mydata[,c(1,2,4)]
table(teste[,3])
cm<-neuralbruno(teste, layers = c(8,6))
confusionMatrix(cm$predicted, cm$original)
#"cat_pH"
teste<-mydata[,c(1,2,5)]
table(teste[,3])
cm<-neuralbruno(teste, layers = 4)
confusionMatrix(cm$predicted, cm$original)
#"Cat.cla"
teste<-mydata[,c(1,2,6)]
table(teste[,3])
cm<-neuralbruno(teste, layers = c(6,5))
confusionMatrix(cm$predicted, cm$original)
#"catestadotrof"
teste<-mydata[,c(1,2,7)]
table(teste[,3])
cm<-neuralbruno(teste, layers = 1)
confusionMatrix(cm$predicted, cm$original)
