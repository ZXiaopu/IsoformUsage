library(tidyverse)
library(caret)

args = commandArgs(trailingOnly=TRUE)
id = as.numeric(args[1])
path = as.character(args[2])

setwd(paste("Isoform/02.Data/", path, sep = ""))
files <- list.files()

geneid <- files[id]
data1 <- read.csv(geneid, sep = "\t", check.names = FALSE,
                  header = TRUE, stringsAsFactors = FALSE)
data1$pop[which(data1$pop %in% c("British", "Finnish", "Tuscan", "Utah"))]<-"EUR"
data1$pop[which(data1$pop == "Yoruba")]<-"YRI"
data1$pop <- factor(data1$pop)

data1<-data1[,-1]

#drop individuals with less than ten reads for now
data1<-data1[which(rowSums(data1[,-c(dim(data1)[2])]) > 10),]
data1[,-c(dim(data1)[2])] <- data1[,-c(dim(data1)[2])]/rowSums(data1[,-c(dim(data1)[2])])

TrainVal <- createDataPartition (y = data1$pop, 
                                 p = 0.80, 
                                 list = FALSE)
training <- data1[TrainVal,]
testing <- data1[-TrainVal,]
ControlParameters <- trainControl(method = "repeatedcv", 
                                  number = 10, repeats = 5,
                                  classProbs = TRUE, summaryFunction = twoClassSummary,
                                  savePredictions = "final", allowParallel = TRUE)

model <- train(pop ~ ., data=training, method="glm",
               preProcess=c("nzv", "center", "scale"),
               trControl=ControlParameters)
pred <- predict(model, testing)
conf <- confusionMatrix(pred, factor(testing$pop))
mat<-t(rbind(as.matrix(conf, what="classes"),as.matrix(conf, what="overall")))
rownames(mat) <- geneid

setwd("Isoform/03.ModelTraining/Baseline")
write.table(mat, paste(path, "baseline", sep = "."),
            col.names = F, row.names = T, quote = F, append = T)
