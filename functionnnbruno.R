# load libs
require(neuralnet)
require(nnet)
require(ggplot2)

# load data

load("/home/bruno/Documentos/DOC/MAIN/100LAGOS/Script&Data/Data100Lakes.RData")
rm(list=ls()[-3]) 

colunas= c(38, 39, 42, 43, 45)
preditores=c(1:3)
(mydata<-dataall[,c(preditores, colunas)])
neuralbruno<-function(mydata, layers=5){
  mydata.<-mydata
  # turnin into factor
  mydata.[,4]<-factor(mydata.[,4])
  # training set
  train<-cbind(mydata.[,-4], class.ind(mydata.[,4]))
  #renamiing
  nlabels<-ncol(train)-3
  
  if(nlabels==3){
    names(train) <- c(names(mydata.)[-4],"l1","l2","l3")
  }
  if(nlabels==4){
    names(train) <- c(names(mydata.)[-4],"l1","l2","l3", "l4")
  }
  if(nlabels==5){
    names(train) <- c(names(mydata.)[-4],"l1","l2","l3", "l4", "l5")
  }
  
  #scaling
  scl <- function(x){ (x - min(x))/(max(x) - min(x)) }
  train[, c(1:3)] <- data.frame(lapply(train[, c(1:3)], scl))
  
  # renaming and formula to neuralnet
  n <- names(train)
  
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
                  lifesign = "minimal", threshold = .02)
  
  #plot(nn)
  # Compute predictions
  pr.nn <- compute(nn, train[, c(1:3)])
  
  # Extract results
  pr.nn_ <- pr.nn$net.result
  head(pr.nn_)
  
  original_values <- max.col(train[, nlabels:ncol(train)])
  pr.nn_2 <- max.col(pr.nn_)
  return(mean(pr.nn_2 == original_values))
}


for( i in 4:ncol(mydata)){
  teste<-mydata[,c(1,2,3,i)]
  print(neuralbruno(teste))
}

