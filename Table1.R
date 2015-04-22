library(RWeka)
library(caret)
library(mlbench)
library(Hmisc)
library(randomForest)
library(kernlab)
library(e1071)
library(cluster) 
library(FactoMineR)
library(randomGLM)
library(corrplot)
library(C50)
library(nnet)
library(e1071)
library(GA)
library(cvTools) 
library(Metrics)
library(MASS)
library(plsdepot)
library(protr)



x = readFASTA("OliFP.Fasta")
x <- x[(sapply(x, protcheck))] # This check unnusal amino acid letters from fasta
length(x)
AAC <- t(sapply(x,extractAAC)) # Amino acid Composition (AAC)
DPC <- t(sapply(x,extractDC)) # Dipeptide Compsition (DPC)
PCP <- t(sapply(x, extractMoreauBroto, props = AAindex$AccNo, nlag = 1L)) # Physiochemical Properties (PCP) from AAindex Database
AAC <- data.frame(AAC, row.names=NULL)
DPC <- data.frame(DPC, row.names=NULL)
PCP <- data.frame(PCP, row.names=NULL)
PCP <- PCP[,colSums(is.na(PCP)) !=nrow(PCP)] # Initial 544 with NA are removed and left with 531 PCP Descriptors

Y = read.csv("Oligomerization.csv", header = TRUE)
class = Y[,2]

Model1 = cbind(AAC,class)
Model2 = cbind(DPC,class)
Model3 = cbind(PCP,class)
Model4 = cbind(AAC,DPC,class)
Model5 = cbind(AAC,PCP,class)
Model6 = cbind(DPC,PCP,class)
Model7 = cbind(AAC,DPC,PCP,class)
## Sampling for 100 Times for model 7
set.seed(4)
m=100 
Resultfull <- matrix(nrow = 4, ncol = m)
Resultcv <- matrix(nrow = 4, ncol = m)
Resultext <- matrix(nrow = 4, ncol = m)

for (i in 1:m){
  
  sample1 <- c(sample(1:194 ,39))
  sample2 <- c(sample(1:203,41))
  Pos <- subset(Model7 , class == 'Oligomer')
  Neg <- subset(Model7 , class == 'Monomer')
  
  train1  <- Pos[-sample1,]
  train2  <- Neg[-sample2,]
  test1 <-   Pos[sample1,]   
  test2 <-   Neg[sample2,]   
  internal <- rbind(train1,train2)
  external <- rbind(test1,test2)
  
  #Training Set
  m2 <- J48(class ~ ., data = internal)
  predcv = table(internal$class, predict(m2,internal[,-ncol(internal)]))
  Resultfull[,i] <- rbind(predcv[1], predcv[3],predcv[2], predcv[4])
  
  # 10-fold Cross Validation
  k=10 #Folds
  id <- sample(1:k,nrow(internal),replace=TRUE)
  list <- 1:k
  Truepos <- data.frame()
  Trueneg <- data.frame()
  Falsepos <- data.frame()
  Falseneg <- data.frame()
  progress.bar <- create_progress_bar("text")
  progress.bar$init(k)
  
  for (h in 1:k){
    train <- subset(internal, id %in% list[-h])
    test <- subset(internal, id %in% c(h))
    m2 <- J48(class ~ ., data = train)
    predcv = table(test$class, predict(m2,test[,-ncol(test)]))
    Truepos <- rbind(Truepos, predcv[4])
    Trueneg <- rbind(Trueneg, predcv[1])
    Falsepos <- rbind(Falsepos, predcv[3])
    Falseneg <- rbind(Falseneg, predcv[2])
    progress.bar$step()
  }
  Resultcv[,i] <- rbind(sum(Trueneg), sum(Falsepos),sum(Falseneg),sum(Truepos))
  
  #External validation
  m2 <- J48(class ~ ., data = internal)
  predcv = table(external$class, predict(m2,external[,-ncol(external)]))
  Resultext[,i] <- rbind(predcv[1], predcv[3],predcv[2], predcv[4])
}

################### Performance report

ACC  <- matrix(nrow = m, ncol = 1)
SENS  <- matrix(nrow = m, ncol = 1)
SPEC  <-matrix(nrow = m, ncol = 1)
MCC <- matrix(nrow = m, ncol = 1)
data = Resultfull
for(i in 1:m){ 
  ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
  SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
  SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
  MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
  MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
  MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
  MCC4  =  sqrt(MCC2)*sqrt(MCC3)
  MCC[i,1]  = MCC1/MCC4}

Perffull = data.frame (mean(ACC),sd(ACC),mean(SENS),sd(SENS),mean(SPEC),sd(SPEC),mean(MCC),sd(MCC))
data = Resultcv
for(i in 1:m){ 
  ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
  SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
  SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
  MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
  MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
  MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
  MCC4	=  sqrt(MCC2)*sqrt(MCC3)
  MCC[i,1]  = MCC1/MCC4}
Perfcv = data.frame (mean(ACC),sd(ACC),mean(SENS),sd(SENS),mean(SPEC),sd(SPEC),mean(MCC),sd(MCC))
data = Resultext
for(i in 1:m){ 
  ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
  SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
  SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
  MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
  MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
  MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
  MCC4	=  sqrt(MCC2)*sqrt(MCC3)
  MCC[i,1]  = MCC1/MCC4}
Perfext = data.frame (mean(ACC),sd(ACC),mean(SENS),sd(SENS),mean(SPEC),sd(SPEC),mean(MCC),sd(MCC))
## record the performacne model 7 as a csv file
write.csv(rbind(Perffull, Perfcv,Perfext), "Performance Model-7.csv", row.names=TRUE, na="")

## Model 6
set.seed(100)
m=100 
Resultfull <- matrix(nrow = 4, ncol = m)
Resultcv <- matrix(nrow = 4, ncol = m)
Resultext <- matrix(nrow = 4, ncol = m)

for (i in 1:m){
  
  sample1 <- c(sample(1:194 ,39))
  sample2 <- c(sample(1:203,41))
  Pos <- subset(Model6 , class == 'Oligomer')
  Neg <- subset(Model6 , class == 'Monomer')
  
  train1  <- Pos[-sample1,]
  train2  <- Neg[-sample2,]
  test1 <-   Pos[sample1,]   
  test2 <-   Neg[sample2,]   
  internal <- rbind(train1,train2)
  external <- rbind(test1,test2)
  
  #Training Set
  m2 <- J48(class ~ ., data = internal)
  predcv = table(internal$class, predict(m2,internal[,-ncol(internal)]))
  Resultfull[,i] <- rbind(predcv[1], predcv[3],predcv[2], predcv[4])
  
  # 10-fold Cross Validation
  k=10 #Folds
  id <- sample(1:k,nrow(internal),replace=TRUE)
  list <- 1:k
  Truepos <- data.frame()
  Trueneg <- data.frame()
  Falsepos <- data.frame()
  Falseneg <- data.frame()
  progress.bar <- create_progress_bar("text")
  progress.bar$init(k)
  
  for (h in 1:k){
    train <- subset(internal, id %in% list[-h])
    test <- subset(internal, id %in% c(h))
    m2 <- J48(class ~ ., data = train)
    predcv = table(test$class, predict(m2,test[,-ncol(test)]))
    Truepos <- rbind(Truepos, predcv[4])
    Trueneg <- rbind(Trueneg, predcv[1])
    Falsepos <- rbind(Falsepos, predcv[3])
    Falseneg <- rbind(Falseneg, predcv[2])
    progress.bar$step()
  }
  Resultcv[,i] <- rbind(sum(Trueneg), sum(Falsepos),sum(Falseneg),sum(Truepos))
  
  #External validation
  m2 <- J48(class ~ ., data = internal)
  predcv = table(external$class, predict(m2,external[,-ncol(external)]))
  Resultext[,i] <- rbind(predcv[1], predcv[3],predcv[2], predcv[4])
}

################### Performance report

ACC  <- matrix(nrow = m, ncol = 1)
SENS  <- matrix(nrow = m, ncol = 1)
SPEC  <-matrix(nrow = m, ncol = 1)
MCC <- matrix(nrow = m, ncol = 1)
data = Resultfull
for(i in 1:m){ 
  ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
  SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
  SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
  MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
  MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
  MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
  MCC4  =  sqrt(MCC2)*sqrt(MCC3)
  MCC[i,1]  = MCC1/MCC4}

Perffull = data.frame (mean(ACC),sd(ACC),mean(SENS),sd(SENS),mean(SPEC),sd(SPEC),mean(MCC),sd(MCC))
data = Resultcv
for(i in 1:m){ 
  ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
  SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
  SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
  MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
  MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
  MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
  MCC4	=  sqrt(MCC2)*sqrt(MCC3)
  MCC[i,1]  = MCC1/MCC4}
Perfcv = data.frame (mean(ACC),sd(ACC),mean(SENS),sd(SENS),mean(SPEC),sd(SPEC),mean(MCC),sd(MCC))
data = Resultext
for(i in 1:m){ 
  ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
  SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
  SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
  MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
  MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
  MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
  MCC4	=  sqrt(MCC2)*sqrt(MCC3)
  MCC[i,1]  = MCC1/MCC4}
Perfext = data.frame (mean(ACC),sd(ACC),mean(SENS),sd(SENS),mean(SPEC),sd(SPEC),mean(MCC),sd(MCC))
## record the performacne model 6 as a csv file
write.csv(rbind(Perffull, Perfcv,Perfext), "Performance Model-6.csv", row.names=TRUE, na="")

## Model 5
set.seed(1)
m=100 
Resultfull <- matrix(nrow = 4, ncol = m)
Resultcv <- matrix(nrow = 4, ncol = m)
Resultext <- matrix(nrow = 4, ncol = m)

for (i in 1:m){
  
  sample1 <- c(sample(1:194 ,39))
  sample2 <- c(sample(1:203,41))
  Pos <- subset(Model5 , class == 'Oligomer')
  Neg <- subset(Model5 , class == 'Monomer')
  
  train1  <- Pos[-sample1,]
  train2  <- Neg[-sample2,]
  test1 <-   Pos[sample1,]   
  test2 <-   Neg[sample2,]   
  internal <- rbind(train1,train2)
  external <- rbind(test1,test2)
  
  #Training Set
  m2 <- J48(class ~ ., data = internal)
  predcv = table(internal$class, predict(m2,internal[,-ncol(internal)]))
  Resultfull[,i] <- rbind(predcv[1], predcv[3],predcv[2], predcv[4])
  
  # 10-fold Cross Validation
  k=10 #Folds
  id <- sample(1:k,nrow(internal),replace=TRUE)
  list <- 1:k
  Truepos <- data.frame()
  Trueneg <- data.frame()
  Falsepos <- data.frame()
  Falseneg <- data.frame()
  progress.bar <- create_progress_bar("text")
  progress.bar$init(k)
  
  for (h in 1:k){
    train <- subset(internal, id %in% list[-h])
    test <- subset(internal, id %in% c(h))
    m2 <- J48(class ~ ., data = train)
    predcv = table(test$class, predict(m2,test[,-ncol(test)]))
    Truepos <- rbind(Truepos, predcv[4])
    Trueneg <- rbind(Trueneg, predcv[1])
    Falsepos <- rbind(Falsepos, predcv[3])
    Falseneg <- rbind(Falseneg, predcv[2])
    progress.bar$step()
  }
  Resultcv[,i] <- rbind(sum(Trueneg), sum(Falsepos),sum(Falseneg),sum(Truepos))
  
  #External validation
  m2 <- J48(class ~ ., data = internal)
  predcv = table(external$class, predict(m2,external[,-ncol(external)]))
  Resultext[,i] <- rbind(predcv[1], predcv[3],predcv[2], predcv[4])
}

################### Performance report

ACC  <- matrix(nrow = m, ncol = 1)
SENS  <- matrix(nrow = m, ncol = 1)
SPEC  <-matrix(nrow = m, ncol = 1)
MCC <- matrix(nrow = m, ncol = 1)
data = Resultfull
for(i in 1:m){ 
  ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
  SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
  SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
  MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
  MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
  MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
  MCC4  =  sqrt(MCC2)*sqrt(MCC3)
  MCC[i,1]  = MCC1/MCC4}

Perffull = data.frame (mean(ACC),sd(ACC),mean(SENS),sd(SENS),mean(SPEC),sd(SPEC),mean(MCC),sd(MCC))
data = Resultcv
for(i in 1:m){ 
  ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
  SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
  SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
  MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
  MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
  MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
  MCC4  =  sqrt(MCC2)*sqrt(MCC3)
  MCC[i,1]  = MCC1/MCC4}
Perfcv = data.frame (mean(ACC),sd(ACC),mean(SENS),sd(SENS),mean(SPEC),sd(SPEC),mean(MCC),sd(MCC))
data = Resultext
for(i in 1:m){ 
  ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
  SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
  SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
  MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
  MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
  MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
  MCC4	=  sqrt(MCC2)*sqrt(MCC3)
  MCC[i,1]  = MCC1/MCC4}
Perfext = data.frame (mean(ACC),sd(ACC),mean(SENS),sd(SENS),mean(SPEC),sd(SPEC),mean(MCC),sd(MCC))
## record the performacne model 5 as a csv file
write.csv(rbind(Perffull, Perfcv,Perfext), "Performance Model-5.csv", row.names=TRUE, na="")

## Model 4
set.seed(5)
m=100 
Resultfull <- matrix(nrow = 4, ncol = m)
Resultcv <- matrix(nrow = 4, ncol = m)
Resultext <- matrix(nrow = 4, ncol = m)

for (i in 1:m){
  
  sample1 <- c(sample(1:194 ,39))
  sample2 <- c(sample(1:203,41))
  Pos <- subset(Model4 , class == 'Oligomer')
  Neg <- subset(Model4 , class == 'Monomer')
  
  train1  <- Pos[-sample1,]
  train2  <- Neg[-sample2,]
  test1 <-   Pos[sample1,]   
  test2 <-   Neg[sample2,]   
  internal <- rbind(train1,train2)
  external <- rbind(test1,test2)
  
  #Training Set
  m2 <- J48(class ~ ., data = internal)
  predcv = table(internal$class, predict(m2,internal[,-ncol(internal)]))
  Resultfull[,i] <- rbind(predcv[1], predcv[3],predcv[2], predcv[4])
  
  # 10-fold Cross Validation
  k=10 #Folds
  id <- sample(1:k,nrow(internal),replace=TRUE)
  list <- 1:k
  Truepos <- data.frame()
  Trueneg <- data.frame()
  Falsepos <- data.frame()
  Falseneg <- data.frame()
  progress.bar <- create_progress_bar("text")
  progress.bar$init(k)
  
  for (h in 1:k){
    train <- subset(internal, id %in% list[-h])
    test <- subset(internal, id %in% c(h))
    m2 <- J48(class ~ ., data = train)
    predcv = table(test$class, predict(m2,test[,-ncol(test)]))
    Truepos <- rbind(Truepos, predcv[4])
    Trueneg <- rbind(Trueneg, predcv[1])
    Falsepos <- rbind(Falsepos, predcv[3])
    Falseneg <- rbind(Falseneg, predcv[2])
    progress.bar$step()
  }
  Resultcv[,i] <- rbind(sum(Trueneg), sum(Falsepos),sum(Falseneg),sum(Truepos))
  
  #External validation
  m2 <- J48(class ~ ., data = internal)
  predcv = table(external$class, predict(m2,external[,-ncol(external)]))
  Resultext[,i] <- rbind(predcv[1], predcv[3],predcv[2], predcv[4])
}

################### Performance report

ACC  <- matrix(nrow = m, ncol = 1)
SENS  <- matrix(nrow = m, ncol = 1)
SPEC  <-matrix(nrow = m, ncol = 1)
MCC <- matrix(nrow = m, ncol = 1)
data = Resultfull
for(i in 1:m){ 
  ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
  SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
  SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
  MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
  MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
  MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
  MCC4  =  sqrt(MCC2)*sqrt(MCC3)
  MCC[i,1]  = MCC1/MCC4}

Perffull = data.frame (mean(ACC),sd(ACC),mean(SENS),sd(SENS),mean(SPEC),sd(SPEC),mean(MCC),sd(MCC))
data = Resultcv
for(i in 1:m){ 
  ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
  SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
  SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
  MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
  MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
  MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
  MCC4  =  sqrt(MCC2)*sqrt(MCC3)
  MCC[i,1]  = MCC1/MCC4}
Perfcv = data.frame (mean(ACC),sd(ACC),mean(SENS),sd(SENS),mean(SPEC),sd(SPEC),mean(MCC),sd(MCC))
data = Resultext
for(i in 1:m){ 
  ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
  SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
  SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
  MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
  MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
  MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
  MCC4	=  sqrt(MCC2)*sqrt(MCC3)
  MCC[i,1]  = MCC1/MCC4}
Perfext = data.frame (mean(ACC),sd(ACC),mean(SENS),sd(SENS),mean(SPEC),sd(SPEC),mean(MCC),sd(MCC))
## record the performacne model 4 as a csv file
write.csv(rbind(Perffull, Perfcv,Perfext), "Performance Model-4.csv", row.names=TRUE, na="")

## Model 3
set.seed(10)
m=100 
Resultfull <- matrix(nrow = 4, ncol = m)
Resultcv <- matrix(nrow = 4, ncol = m)
Resultext <- matrix(nrow = 4, ncol = m)

for (i in 1:m){
  
  sample1 <- c(sample(1:194 ,39))
  sample2 <- c(sample(1:203,41))
  Pos <- subset(Model3 , class == 'Oligomer')
  Neg <- subset(Model3 , class == 'Monomer')
  
  train1  <- Pos[-sample1,]
  train2  <- Neg[-sample2,]
  test1 <-   Pos[sample1,]   
  test2 <-   Neg[sample2,]   
  internal <- rbind(train1,train2)
  external <- rbind(test1,test2)
  
  #Training Set
  m2 <- J48(class ~ ., data = internal)
  predcv = table(internal$class, predict(m2,internal[,-ncol(internal)]))
  Resultfull[,i] <- rbind(predcv[1], predcv[3],predcv[2], predcv[4])
  
  # 10-fold Cross Validation
  k=10 #Folds
  id <- sample(1:k,nrow(internal),replace=TRUE)
  list <- 1:k
  Truepos <- data.frame()
  Trueneg <- data.frame()
  Falsepos <- data.frame()
  Falseneg <- data.frame()
  progress.bar <- create_progress_bar("text")
  progress.bar$init(k)
  
  for (h in 1:k){
    train <- subset(internal, id %in% list[-h])
    test <- subset(internal, id %in% c(h))
    m2 <- J48(class ~ ., data = train)
    predcv = table(test$class, predict(m2,test[,-ncol(test)]))
    Truepos <- rbind(Truepos, predcv[4])
    Trueneg <- rbind(Trueneg, predcv[1])
    Falsepos <- rbind(Falsepos, predcv[3])
    Falseneg <- rbind(Falseneg, predcv[2])
    progress.bar$step()
  }
  Resultcv[,i] <- rbind(sum(Trueneg), sum(Falsepos),sum(Falseneg),sum(Truepos))
  
  #External validation
  m2 <- J48(class ~ ., data = internal)
  predcv = table(external$class, predict(m2,external[,-ncol(external)]))
  Resultext[,i] <- rbind(predcv[1], predcv[3],predcv[2], predcv[4])
}

################### Performance report

ACC  <- matrix(nrow = m, ncol = 1)
SENS  <- matrix(nrow = m, ncol = 1)
SPEC  <-matrix(nrow = m, ncol = 1)
MCC <- matrix(nrow = m, ncol = 1)
data = Resultfull
for(i in 1:m){ 
  ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
  SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
  SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
  MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
  MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
  MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
  MCC4  =  sqrt(MCC2)*sqrt(MCC3)
  MCC[i,1]  = MCC1/MCC4}

Perffull = data.frame (mean(ACC),sd(ACC),mean(SENS),sd(SENS),mean(SPEC),sd(SPEC),mean(MCC),sd(MCC))
data = Resultcv
for(i in 1:m){ 
  ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
  SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
  SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
  MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
  MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
  MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
  MCC4  =  sqrt(MCC2)*sqrt(MCC3)
  MCC[i,1]  = MCC1/MCC4}
Perfcv = data.frame (mean(ACC),sd(ACC),mean(SENS),sd(SENS),mean(SPEC),sd(SPEC),mean(MCC),sd(MCC))
data = Resultext
for(i in 1:m){ 
  ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
  SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
  SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
  MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
  MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
  MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
  MCC4	=  sqrt(MCC2)*sqrt(MCC3)
  MCC[i,1]  = MCC1/MCC4}
Perfext = data.frame (mean(ACC),sd(ACC),mean(SENS),sd(SENS),mean(SPEC),sd(SPEC),mean(MCC),sd(MCC))
## record the performacne model 3 as a csv file
write.csv(rbind(Perffull, Perfcv,Perfext), "Performance Model-3.csv", row.names=TRUE, na="")

## Model 2
set.seed(24)
m=100 
Resultfull <- matrix(nrow = 4, ncol = m)
Resultcv <- matrix(nrow = 4, ncol = m)
Resultext <- matrix(nrow = 4, ncol = m)

for (i in 1:m){
  
  sample1 <- c(sample(1:194 ,39))
  sample2 <- c(sample(1:203,41))
  Pos <- subset(Model2 , class == 'Oligomer')
  Neg <- subset(Model2 , class == 'Monomer')
  
  train1  <- Pos[-sample1,]
  train2  <- Neg[-sample2,]
  test1 <-   Pos[sample1,]   
  test2 <-   Neg[sample2,]   
  internal <- rbind(train1,train2)
  external <- rbind(test1,test2)
  
  #Training Set
  m2 <- J48(class ~ ., data = internal)
  predcv = table(internal$class, predict(m2,internal[,-ncol(internal)]))
  Resultfull[,i] <- rbind(predcv[1], predcv[3],predcv[2], predcv[4])
  
  # 10-fold Cross Validation
  k=10 #Folds
  id <- sample(1:k,nrow(internal),replace=TRUE)
  list <- 1:k
  Truepos <- data.frame()
  Trueneg <- data.frame()
  Falsepos <- data.frame()
  Falseneg <- data.frame()
  progress.bar <- create_progress_bar("text")
  progress.bar$init(k)
  
  for (h in 1:k){
    train <- subset(internal, id %in% list[-h])
    test <- subset(internal, id %in% c(h))
    m2 <- J48(class ~ ., data = train)
    predcv = table(test$class, predict(m2,test[,-ncol(test)]))
    Truepos <- rbind(Truepos, predcv[4])
    Trueneg <- rbind(Trueneg, predcv[1])
    Falsepos <- rbind(Falsepos, predcv[3])
    Falseneg <- rbind(Falseneg, predcv[2])
    progress.bar$step()
  }
  Resultcv[,i] <- rbind(sum(Trueneg), sum(Falsepos),sum(Falseneg),sum(Truepos))
  
  #External validation
  m2 <- J48(class ~ ., data = internal)
  predcv = table(external$class, predict(m2,external[,-ncol(external)]))
  Resultext[,i] <- rbind(predcv[1], predcv[3],predcv[2], predcv[4])
}

################### Performance report

ACC  <- matrix(nrow = m, ncol = 1)
SENS  <- matrix(nrow = m, ncol = 1)
SPEC  <-matrix(nrow = m, ncol = 1)
MCC <- matrix(nrow = m, ncol = 1)
data = Resultfull
for(i in 1:m){ 
  ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
  SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
  SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
  MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
  MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
  MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
  MCC4  =  sqrt(MCC2)*sqrt(MCC3)
  MCC[i,1]  = MCC1/MCC4}

Perffull = data.frame (mean(ACC),sd(ACC),mean(SENS),sd(SENS),mean(SPEC),sd(SPEC),mean(MCC),sd(MCC))
data = Resultcv
for(i in 1:m){ 
  ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
  SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
  SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
  MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
  MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
  MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
  MCC4  =  sqrt(MCC2)*sqrt(MCC3)
  MCC[i,1]  = MCC1/MCC4}
Perfcv = data.frame (mean(ACC),sd(ACC),mean(SENS),sd(SENS),mean(SPEC),sd(SPEC),mean(MCC),sd(MCC))
data = Resultext
for(i in 1:m){ 
  ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
  SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
  SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
  MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
  MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
  MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
  MCC4	=  sqrt(MCC2)*sqrt(MCC3)
  MCC[i,1]  = MCC1/MCC4}
Perfext = data.frame (mean(ACC),sd(ACC),mean(SENS),sd(SENS),mean(SPEC),sd(SPEC),mean(MCC),sd(MCC))
## record the performacne model 2 as a csv file
write.csv(rbind(Perffull, Perfcv,Perfext), "Performance Model-2.csv", row.names=TRUE, na="")

## Model 1
set.seed(99)
m=100 
Resultfull <- matrix(nrow = 4, ncol = m)
Resultcv <- matrix(nrow = 4, ncol = m)
Resultext <- matrix(nrow = 4, ncol = m)

for (i in 1:m){
  
  sample1 <- c(sample(1:194 ,39))
  sample2 <- c(sample(1:203,41))
  Pos <- subset(Model1 , class == 'Oligomer')
  Neg <- subset(Model1 , class == 'Monomer')
  
  train1  <- Pos[-sample1,]
  train2  <- Neg[-sample2,]
  test1 <-   Pos[sample1,]   
  test2 <-   Neg[sample2,]   
  internal <- rbind(train1,train2)
  external <- rbind(test1,test2)
  
  #Training Set
  m2 <- J48(class ~ ., data = internal)
  predcv = table(internal$class, predict(m2,internal[,-ncol(internal)]))
  Resultfull[,i] <- rbind(predcv[1], predcv[3],predcv[2], predcv[4])
  
  # 10-fold Cross Validation
  k=10 #Folds
  id <- sample(1:k,nrow(internal),replace=TRUE)
  list <- 1:k
  Truepos <- data.frame()
  Trueneg <- data.frame()
  Falsepos <- data.frame()
  Falseneg <- data.frame()
  progress.bar <- create_progress_bar("text")
  progress.bar$init(k)
  
  for (h in 1:k){
    train <- subset(internal, id %in% list[-h])
    test <- subset(internal, id %in% c(h))
    m2 <- J48(class ~ ., data = train)
    predcv = table(test$class, predict(m2,test[,-ncol(test)]))
    Truepos <- rbind(Truepos, predcv[4])
    Trueneg <- rbind(Trueneg, predcv[1])
    Falsepos <- rbind(Falsepos, predcv[3])
    Falseneg <- rbind(Falseneg, predcv[2])
    progress.bar$step()
  }
  Resultcv[,i] <- rbind(sum(Trueneg), sum(Falsepos),sum(Falseneg),sum(Truepos))
  
  #External validation
  m2 <- J48(class ~ ., data = internal)
  predcv = table(external$class, predict(m2,external[,-ncol(external)]))
  Resultext[,i] <- rbind(predcv[1], predcv[3],predcv[2], predcv[4])
}

################### Performance report

ACC  <- matrix(nrow = m, ncol = 1)
SENS  <- matrix(nrow = m, ncol = 1)
SPEC  <-matrix(nrow = m, ncol = 1)
MCC <- matrix(nrow = m, ncol = 1)
data = Resultfull
for(i in 1:m){ 
  ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
  SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
  SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
  MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
  MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
  MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
  MCC4  =  sqrt(MCC2)*sqrt(MCC3)
  MCC[i,1]  = MCC1/MCC4}

Perffull = data.frame (mean(ACC),sd(ACC),mean(SENS),sd(SENS),mean(SPEC),sd(SPEC),mean(MCC),sd(MCC))
data = Resultcv
for(i in 1:m){ 
  ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
  SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
  SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
  MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
  MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
  MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
  MCC4  =  sqrt(MCC2)*sqrt(MCC3)
  MCC[i,1]  = MCC1/MCC4}
Perfcv = data.frame (mean(ACC),sd(ACC),mean(SENS),sd(SENS),mean(SPEC),sd(SPEC),mean(MCC),sd(MCC))
data = Resultext
for(i in 1:m){ 
  ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
  SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
  SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
  MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
  MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
  MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
  MCC4	=  sqrt(MCC2)*sqrt(MCC3)
  MCC[i,1]  = MCC1/MCC4}
Perfext = data.frame (mean(ACC),sd(ACC),mean(SENS),sd(SENS),mean(SPEC),sd(SPEC),mean(MCC),sd(MCC))
## record the performacne model 1 as a csv file
write.csv(rbind(Perffull, Perfcv,Perfext), "Performance Model-1.csv", row.names=TRUE, na="")