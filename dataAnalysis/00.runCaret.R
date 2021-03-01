library(tidyverse)
library(caret)
library(xgboost)
library(caretEnsemble)

args = commandArgs(trailingOnly=TRUE)
id = as.numeric(args[1])
path = as.character(args[3])
setwd(paste("Isoform/02.Data", path, sep="/"))
files <- list.files()

print(id)
print(path)

geneid <- files[id]
g <- gsub(paste(".",path,sep=""), "", geneid)
data1 <- read.csv(geneid, sep = "\t", check.names = FALSE,
                  header = TRUE, stringsAsFactors = FALSE)
data1$pop[which(data1$pop %in% c("British", "Finnish", "Tuscan", "Utah"))]<-"EUR"
data1$pop[which(data1$pop == "Yoruba")]<-"YRI"
data1$pop <- factor(data1$pop)

data1<-data1[which(rowSums(data1[,-c(dim(data1)[2])]) > 10),]
data1[,-c(dim(data1)[2])] <- data1[,-c(dim(data1)[2])]/rowSums(data1[,-c(dim(data1)[2])])

len <- read.csv(paste("Isoform/01.OwnReference/02.MappingModels/win3.flank1500/", ".len", sep=g), sep="\t", header=F)
len$length <- (len$V3-len$V2)+1

data1[,-c(dim(data1)[2])] <-data.frame(t(t(data1[,-c(dim(data1)[2])])/len$length))

if(length(table(factor(data1$pop)))==1)
  {print("skip")}else{
  
    if(length(nearZeroVar(data1))>0 & ncol(data1[,-nearZeroVar(data1)])<3) 
    {print("skip")} else{
    
    TrainVal <- createDataPartition (y = data1$pop, 
                                     p = 0.80, 
                                     list = FALSE)
    training <- data1[TrainVal,]
    testing <- data1[-TrainVal,]
    #if dont want to permute population column set second argument to 0
    #otherwise this specifies the permutation number
    if(args[2] > 0)
    {
      set.seed(args[2])
      training$pop<-sample(training$pop, replace = FALSE)
    }
    
    ControlParamteres <- trainControl(method = "repeatedcv", 
                                      number = 10, repeats = 5,
                                      classProbs = TRUE, summaryFunction = twoClassSummary,
                                      savePredictions = "final", allowParallel = TRUE)
    xgbGrid <- expand.grid(eta=c(0.01,0.05,0.1),
                           max_depth=c(2,4,6,8),
                           nrounds=c(50,100,150),
                           min_child_weight=1,
                           subsample=c(0.5,0.75,1.0),
                           colsample_bytree=c(0.4,0.6,0.8),
                           gamma=0)
    models <- caretList(pop ~ ., data = training, trControl = ControlParamteres, metric = "ROC",
                        tuneList = list(logit = caretModelSpec(method = "glm", family = "binomial"),
                                        elasticnet = caretModelSpec(method = "glmnet",
                                                                    tuneGrid = expand.grid(alpha = 0:1,
                                                                                           lambda = seq(0.0001, 1, length = 20))),
                                        rf = caretModelSpec(method = "ranger"),
                                        xgbmodel = caretModelSpec(method = "xgbTree",tuneGrid=xgbGrid, nthread=1)),
                        preProcess = c("nzv", "center", "scale"))
    models.preds <- lapply(models, predict, newdata = testing) #add type = "prob" for class probabilities
    models.preds <- data.frame(models.preds)
    
    results<-list()
    for(i in 1:4)
    {
      conf<-confusionMatrix(models.preds[,i], as.factor(testing$pop))
      mat<-t(rbind(features=dim(data1)[2],as.matrix(conf, what="classes"),as.matrix(conf, what="overall")))
      results[[paste(g, i, sep=" ")]]<-mat
    }
    res <- data.frame(results %>% map_df(as_tibble))
    rownames(res) <-paste(g,rownames(res), sep=" ")
    colnames(res)<-colnames(results[[1]])
    setwd(paste("Isoform/03.ModelTraining/01.metrics.20201028", paste("metrics",args[2],sep=""), sep="/"))
    write.table(res, paste(g,"_metrics.txt", sep=""))
  }}


