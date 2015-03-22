#!/usr/bin/env Rscript
library(RCurl)
set.seed(3)
x <- getURL("https://raw.githubusercontent.com/Rnewbie/OliFP/master/OliFP.csv")
OliFP <- read.csv(text=x, header = TRUE)
data <- OliFP
PCP <- data[, 2:532]
DPC <- data[, 533:932]
AAC <- data[, 933:952]
AAC_DPC <- cbind(AAC, DPC)
AAC_PCP <- cbind(AAC, PCP)
DPC_PCP <- cbind(DPC, PCP)
ALL <- data[, 5:952]
Oligomerization <- data$Oligomerization
set.seed(1)
## Bind the Oligomerization
PCP <- cbind(Oligomerization, PCP)
DPC <- cbind(Oligomerization, DPC)
AAC <- cbind(Oligomerization, AAC)
AAC_DPC <- cbind(Oligomerization, AAC_DPC)
AAC_PCP <- cbind(Oligomerization, AAC_PCP)
DPC_PCP <- cbind(Oligomerization, DPC_PCP)
ALL <- cbind(Oligomerization, ALL)
## Seperate the Oligomer and Monomer
set.seed(1)
library(caret)
PCP_Monomer <- subset(PCP, Oligomerization == "Monomer")
PCP_Oligomer <- subset(PCP, Oligomerization == "Oligomer")
DPC_Monomer <- subset(DPC, Oligomerization == "Monomer")
DPC_Oligomer <- subset(DPC, Oligomerization == "Oligomer")
AAC_Monomer <- subset(AAC, Oligomerization == "Monomer")
AAC_Oligomer <- subset(AAC, Oligomerization == "Oligomer")
AAC_DPC_Monomer <- subset(AAC_DPC, Oligomerization == "Monomer")
AAC_DPC_Oligomer <- subset(AAC_DPC, Oligomerization == "Oligomer")
AAC_PCP_Monomer <- subset(AAC_PCP, Oligomerization == "Monomer")
AAC_PCP_Oligomer <- subset(AAC_PCP, Oligomerization == "Oligomer")
DPC_PCP_Monomer <- subset(DPC_PCP, Oligomerization == "Monomer")
DPC_PCP_Oligomer <- subset(DPC_PCP, Oligomerization == "Oligomer")
ALL_Monomer <- subset(ALL, Oligomerization == "Monomer")
ALL_Oligomer <- subset(ALL, Oligomerization == "Oligomer")
## Do Kennard Stone Data Partationing
library(prospectr)
set.seed(6)
x <- data.frame(PCP_Monomer)
sel <- kenStone(x[-1], k = 150, metric = "mahal", pc=2)
trainPCP_Monomer <- x[sel$model, ]
testPCP_Monomer <- x[sel$test, ]
x <- data.frame(PCP_Oligomer)
sel <- kenStone(x[-1], k = 149, metric = "mahal", pc=2)
trainPCP_Oligomer <- x[sel$model, ]
testPCP_Oligomer <- x[sel$test, ]

x <- data.frame(DPC_Monomer)
sel <- kenStone(x[-1], k = 150, metric = "mahal", pc=2)
trainDPC_Monomer <- x[sel$model, ]
testDPC_Monomer <- x[sel$test, ]
x <- data.frame(DPC_Oligomer)
sel <- kenStone(x[-1], k = 149, metric = "mahal", pc=2)
trainDPC_Oligomer <- x[sel$model, ]
testDPC_Oligomer <- x[sel$test, ]

x <- data.frame(AAC_Monomer)
sel <- kenStone(x[-1], k = 150, metric = "mahal", pc=2)
trainAAC_Monomer <- x[sel$model, ]
testAAC_Monomer <- x[sel$test, ]
x <- data.frame(AAC_Oligomer)
sel <- kenStone(x[-1], k = 149, metric = "mahal", pc=2)
trainAAC_Oligomer <- x[sel$model, ]
testAAC_Oligomer <- x[sel$test, ]

x <- data.frame(AAC_DPC_Monomer)
sel <- kenStone(x[-1], k = 150, metric = "mahal", pc=2)
trainAAC_DPC_Monomer <- x[sel$model, ]
testAAC_DPC_Monomer <- x[sel$test, ]
x <- data.frame(AAC_DPC_Oligomer)
sel <- kenStone(x[-1], k = 149, metric = "mahal", pc=2)
trainAAC_DPC_Oligomer <- x[sel$model, ]
testAAC_DPC_Oligomer <- x[sel$test, ]

x <- data.frame(AAC_PCP_Monomer)
sel <- kenStone(x[-1], k = 150, metric = "mahal", pc=2)
trainAAC_PCP_Monomer <- x[sel$model, ]
testAAC_PCP_Monomer <- x[sel$test, ]
x <- data.frame(AAC_PCP_Oligomer)
sel <- kenStone(x[-1], k = 149, metric = "mahal", pc=2)
trainAAC_PCP_Oligomer <- x[sel$model, ]
testAAC_PCP_Oligomer <- x[sel$test, ]

x <- data.frame(DPC_PCP_Monomer)
sel <- kenStone(x[-1], k = 150, metric = "mahal", pc=2)
trainDPC_PCP_Monomer <- x[sel$model, ]
testDPC_PCP_Monomer <- x[sel$test, ]
x <- data.frame(DPC_PCP_Oligomer)
sel <- kenStone(x[-1], k = 149, metric = "mahal", pc=2)
trainDPC_PCP_Oligomer <- x[sel$model, ]
testDPC_PCP_Oligomer <- x[sel$test, ]

x <- data.frame(ALL_Monomer)
sel <- kenStone(x[-1], k = 150, metric = "mahal", pc=2)
trainALL_Monomer <- x[sel$model, ]
testALL_Monomer <- x[sel$test, ]
x <- data.frame(ALL_Oligomer)
sel <- kenStone(x[-1], k = 149, metric = "mahal", pc=2)
trainALL_Oligomer <- x[sel$model, ]
testALL_Oligomer <- x[sel$test, ]
## bind Oligomer and Monomer or training set and testing set
set.seed(400)
PCP_Train <- rbind(trainPCP_Monomer, trainPCP_Oligomer)
PCP_Test <- rbind(testPCP_Monomer, testPCP_Oligomer)
DPC_Train <- rbind(trainDPC_Monomer, trainDPC_Oligomer)
DPC_Test <- rbind(testDPC_Monomer, testDPC_Oligomer)
AAC_Train <- rbind(trainAAC_Monomer, trainAAC_Oligomer)
AAC_Test <- rbind(testAAC_Monomer, testAAC_Oligomer)
AAC_DPC_Train <- rbind(trainAAC_DPC_Monomer, trainAAC_DPC_Oligomer)
AAC_DPC_Test <- rbind(testAAC_DPC_Monomer, testAAC_DPC_Oligomer)
AAC_PCP_Train <- rbind(trainAAC_PCP_Monomer, trainAAC_PCP_Oligomer)
AAC_PCP_Test <- rbind(testAAC_PCP_Monomer, testAAC_PCP_Oligomer)
DPC_PCP_Train <- rbind(trainDPC_PCP_Monomer, trainDPC_PCP_Oligomer)
DPC_PCP_Test <- rbind(testDPC_PCP_Monomer, testDPC_PCP_Oligomer)
ALL_Train <- rbind(trainALL_Monomer, trainALL_Oligomer)
ALL_Test <- rbind(testALL_Monomer, testALL_Oligomer)
## Building Predictive Model for PCP, training set, 10-fold and testing set
library(RWeka)
set.seed(1)
PCPfit <- J48(Oligomerization~., data = PCP_Train)
summary(PCPfit)
eval10fold <- evaluate_Weka_classifier(PCPfit,
                                       numFolds=10,
                                       complexity = FALSE,
                                       seed=1,
                                       class=TRUE)
eval10fold
external <- evaluate_Weka_classifier(PCPfit,
                                     newdata = PCP_Test,
                                       numFolds=10,
                                       complexity = FALSE,
                                       seed=1,
                                       class=TRUE)
external
## Building Predictive Model for DPC, training set, 10-fold and testing set
set.seed(2)
DPCfit <- J48(Oligomerization~., data = DPC_Train)
summary(DPCfit)
eval10fold <- evaluate_Weka_classifier(DPCfit,
                                       numFolds=10,
                                       complexity = FALSE,
                                       seed=1,
                                       class=TRUE)
eval10fold
external <- evaluate_Weka_classifier(DPCfit,
                                     newdata = DPC_Test,
                                     numFolds=10,
                                     complexity = FALSE,
                                     seed=1,
                                     class=TRUE)
external

## Building Predictive Model for DPC, training set, 10-fold and testing set
set.seed(3)
AACfit <- J48(Oligomerization~., data = AAC_Train)
summary(AACfit)
eval10fold <- evaluate_Weka_classifier(AACfit,
                                       numFolds=10,
                                       complexity = FALSE,
                                       seed=1,
                                       class=TRUE)
eval10fold
external <- evaluate_Weka_classifier(AACfit,
                                     newdata = AAC_Test,
                                     numFolds=10,
                                     complexity = FALSE,
                                     seed=1,
                                     class=TRUE)
external

## Building Predictive Model for AAC_DPC, training set, 10-fold and testing set
set.seed(4)
AAC_DPCfit <- J48(Oligomerization~., data = AAC_DPC_Train)
summary(AAC_DPCfit)
eval10fold <- evaluate_Weka_classifier(AAC_DPCfit,
                                       numFolds=10,
                                       complexity = FALSE,
                                       seed=1,
                                       class=TRUE)
eval10fold
external <- evaluate_Weka_classifier(AAC_DPCfit,
                                     newdata = AAC_DPC_Test,
                                     numFolds=10,
                                     complexity = FALSE,
                                     seed=1,
                                     class=TRUE)
external

## Building Predictive Model for AAC_PCP, training set, 10-fold and testing set
set.seed(5)
AAC_PCPfit <- J48(Oligomerization~., data = AAC_PCP_Train)
summary(AAC_PCPfit)
eval10fold <- evaluate_Weka_classifier(AAC_PCPfit,
                                       numFolds=10,
                                       complexity = FALSE,
                                       seed=1,
                                       class=TRUE)
eval10fold
external <- evaluate_Weka_classifier(AAC_PCPfit,
                                     newdata = AAC_PCP_Test,
                                     numFolds=10,
                                     complexity = FALSE,
                                     seed=1,
                                     class=TRUE)
external
## Building Predictive Model for DPC_PCP, training set, 10-fold and testing set
set.seed(6)
DPC_PCPfit <- J48(Oligomerization~., data = DPC_PCP_Train)
summary(DPC_PCPfit)
eval10fold <- evaluate_Weka_classifier(DPC_PCPfit,
                                       numFolds=10,
                                       complexity = FALSE,
                                       seed=1,
                                       class=TRUE)
eval10fold
external <- evaluate_Weka_classifier(DPC_PCPfit,
                                     newdata = DPC_PCP_Test,
                                     numFolds=10,
                                     complexity = FALSE,
                                     seed=1,
                                     class=TRUE)
external

## Building Predictive Model for ALL, training set, 10-fold and testing set
set.seed(7)
ALLfit <- J48(Oligomerization~., data = ALL_Train)
summary(ALLfit)
eval10fold <- evaluate_Weka_classifier(ALLfit,
                                       numFolds=10,
                                       complexity = FALSE,
                                       seed=1,
                                       class=TRUE)
eval10fold
external <- evaluate_Weka_classifier(ALLfit,
                                     newdata = ALL_Test,
                                     numFolds=10,
                                     complexity = FALSE,
                                     seed=1,
                                     class=TRUE)
external



#### this is for automation
## automation 
AAC <- read.csv("AAC.csv", header = TRUE)
DPC <- read.csv("DPC.csv", header = TRUE)
PCP <- read.csv("PCP.csv", header = TRUE)
DPC <- data.frame(DPC)
PCP <- data.frame(PCP)
DPC_PCP <- cbind(DPC, PCP)
Oligomerization <- read.csv("Oligomerization.csv", header = TRUE)
Oligomerization <- Oligomerization$Oligomerization
x <- list(DPC_PCP = DPC_PCP)
Logistic <- lapply(x, function(x){
  data <- cbind(Oligomerization, x)
  Monomer <- subset(data, Oligomerization == "Monomer")
  Oligomer <- subset(data, Oligomerization == "Oligomer")
  sel <- kenStone(Monomer[-1], k = 162, metric = "mahal", pc=2)
  train_Monomer <- Monomer[sel$model, ]
  test_Monomer <- Monomer[sel$test, ]
  sel <- kenStone(Oligomer[-1], k = 155, metric = "mahal", pc=2)
  train_Oligomer <- Oligomer[sel$model, ]
  test_Oligomer <- Oligomer[sel$test, ]
  Train <- rbind(train_Monomer, train_Oligomer)
  Test <- rbind(test_Monomer, test_Oligomer)
  #x <- LMT(Oligomerization~., data = Train)
  #cv <- evaluate_Weka_classifier(x,
  #                              numFolds=10,
  #                               complexity = FALSE,
  #                              seed=1,
  #                               class=TRUE)
  #external <- evaluate_Weka_classifier(x,
  #                                    newdata = Test,
  #                                   numFolds=10,
  #                                  complexity = FALSE,
  #                                   seed=1,
  #                                   class=TRUE)
  Model <- summary(x)
  results <- list(Train = Train, Test = Test)
  return(results)
})

Train <- Logistic$DPC_PCP$Train
Test <- Logistic$DPC_PCP$Test
folds <- createFolds(Train$Oligomerization, k = 10)
cv <- lapply(folds, function(x) {
  train <- Train[x, ]
  test <- Train[-x, ]
  model <- randomForest(Oligomerization~., data = train, netree= 500, mtry = 30)
  pred <- predict(model, test)
  actual <- test$Oligomerization
  matrix <- table(pred, actual)
  unlist <- unlist(matrix)
  return(unlist)
})
unlist <- unlist(cv)
df <- data.frame(cv)
x <- unlist(cv$Fold01)
y <- unlist(cv$Fold02)
x <- data.frame(x)
y <- data.frame(y)
results <- cbind(x$Freq, y$Freq)
data <- data.frame(results)
m = ncol(data)
ACC  <- matrix(nrow = m, ncol = 1)
SENS  <- matrix(nrow = m, ncol = 1)
SPEC  <-matrix(nrow = m, ncol = 1)
MCC <- matrix(nrow = m, ncol = 1)

for(i in 1:m){ 
  ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
  SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
  SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
  MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
  MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
  MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
  MCC4  =  sqrt(MCC2)*sqrt(MCC3)
  
  
  MCC[i,1]  = MCC1/MCC4
}


Perf = data.frame (ACC,SENS,SPEC,MCC)
Perf
kable(Perf,
      align = 'c',
      format = "pandoc",
      caption = "Studip Table")
#extra



write.csv(top10AAC, file = "top20AAC.csv")
write.csv(top10DPC, file = "top20DPC.csv")
write.csv(top10PCP, file = "top20PCP.csv")
top10AAC <- read.csv("top20AAC.csv", header = TRUE)
top10DPC <- read.csv("top30DPC.csv", header = TRUE)
top10PCP <- read.csv("top30PCP.csv", header = TRUE)
a <- data.frame(top10AAC)
b <- data.frame(top10DPC)
c <- data.frame(top10PCP)
a$AAC <- factor(a$AAC, levels =a[order(a$Overall), "AAC"])
b$DPC <- factor(b$DPC, levels =b[order(b$Overall), "DPC"])
c$PCP <- factor(c$PCP, levels=c[order(c$Overall), "PCP"])
x <- ggplot(a, aes(x= Overall, y = AAC)) +
  geom_point(size=5, colour = "red") + coord_fixed(14.7) +
  theme_bw() + ggtitle("AAC") + xlab("Feature Usage") + ylab("") +
  scale_x_continuous(breaks = round(seq(min(0), max(100), by = 25),1), limits= c(0, 100))
y <- ggplot(b, aes(x=Overall, y=DPC)) +
  geom_point(size=5, colour = "blue") + coord_fixed(10) +
  theme_bw() + ggtitle("DPC") + xlab("Feature Usage") + ylab("") +
  scale_x_continuous(breaks = round(seq(min(0), max(100), by = 25),1), limits = c(0, 100))
z <- ggplot(c, aes(x=Overall, y=PCP)) +
  geom_point(size=5, colour= "green") + coord_fixed(11) +
  theme_bw() + ggtitle("PCP") + xlab("Feature Usage") + ylab("") +
  scale_x_continuous(breaks = round(seq(min(0), max(100), by = 25),1), limits = c(0, 100))
grid.arrange(x, y, z, ncol = 3)



###
###C5imp(ruleModel, metric = "splits")
