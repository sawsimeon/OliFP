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

