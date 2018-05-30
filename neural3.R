# Script utilizado para modelar as nn's

# load libs
require(neuralnet)
require(nnet)
require(ggplot2)

# load data
rm(list=ls())
load("/home/bruno/Documentos/DOC/MAIN/100LAGOS/Script&Data/Data100Lakes.RData")
rm(list=ls()[-3]) 

colunas= c(38, 39, 42, 43, 45)
preditores=c(1,3)
# 1 e 3
#(mydata<-dataall[,c(preditores, colunas)])
(mydata<-dataall[,c(preditores, colunas)])


neuralbruno<-function(mydata, layers=5){
  mydata.<-mydata
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
  
  plot(nn)
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
  
  # Set seed for reproducibility purposes
  set.seed(500)
  # 10 fold cross validation
  k <- 10
  # Results from cv
  outs <- NULL
  # Train test split proportions
  proportion <- 0.95 # Set to 0.995 for LOOCV
  
  # Crossvalidate, go!
  for(i in 1:k)
  {
    index <- sample(1:nrow(train), round(proportion*nrow(train)))
    train_cv <- train[index, ]
    test_cv <- train[-index, ]
    nn_cv <- neuralnet(f,
                       data = train_cv,
                       hidden = layers,
                       act.fct = "logistic",
                       linear.output = FALSE, threshold = .05)
    
    # Compute predictions
    pr.nn <- compute(nn_cv, test_cv[,1:length(preditores)])
    # Extract results
    pr.nn_ <- pr.nn$net.result
    # Accuracy (test set)
    original_values <- max.col(test_cv[, c((ncol(train)-length(preditores)):ncol(train))])
    pr.nn_2 <- max.col(pr.nn_)
    outs[i] <- mean(pr.nn_2 == original_values)
  }
  
  
  ######################################################
  print(colnames(mydata.))
  print(paste("Run 1: ", media))
  print(paste("Cross: ", mean(outs)))
  
}

for( i in 1:length(colunas)){
  teste<-dataall[,c(preditores, colunas[i])]
  #teste<-mydata[,c(preditores,i)]
  neuralbruno(teste)
}



