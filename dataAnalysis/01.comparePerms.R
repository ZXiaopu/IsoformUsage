library(tidyverse)
library(reshape2)
library(ggpubr)
data_path <- "Isoform/Permutation/"   # path to the data
args = commandArgs(trailingOnly=TRUE)

#######################################################################
#this bit reads the output metrics files from the model fits
readFiles<-function(files, perm)
{
  
  data <- data_frame(filename = files) %>% # create a data frame
    # holding the file names
    mutate(file_contents = map(filename,          # read files into
                               ~ read_delim(file.path(data_path, .), " ", skip=1, col_names = c("gene", "features", "Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value", "Precision", "Recall", "F1", "Prevalence", "Detection Rate", "Detection Prevalence", "Balanced Accuracy", "Accuracy", "Kappa", "AccuracyLower", "AccuracyUpper", "AccuracyNull", "AccuracyPValue", "McnemarPValue"))) # a new data column
    )  
  
  data2<-unnest(data)
  data2<-data2 %>% separate(gene, c("gene", "model"), sep=" ")
  data2$perm<-perm
  save(data2, file=paste("Isoform/Permutation/",perm,".RData", sep=""))
  return(data2)
}

#currently just reads the real model results and results from just one permutation
#however ideally we would run at least ten permutations
files1 <- dir(data_path, pattern = "*[-|+]_metrics.txt") # get file names
datR<-readFiles(files1, "REAL")
files1 <- dir(data_path, pattern = "*1_metrics.txt") # get file names
dat1<-readFiles(files1, number)
files2 <- dir(data_path, pattern = "*2_metrics.txt") # get file names
dat2<-readFiles(files2, 2)
files3 <- dir(data_path, pattern = "*3_metrics.txt") # get file names
dat3<-readFiles(files3, 3)
files4 <- dir(data_path, pattern = "*4_metrics.txt") # get file names
dat4<-readFiles(files4, 4)
files5 <- dir(data_path, pattern = "*5_metrics.txt") # get file names
dat5<-readFiles(files5, 5)
files6 <- dir(data_path, pattern = "*6_metrics.txt") # get file names
dat6<-readFiles(files6, 6)
files7 <- dir(data_path, pattern = "*7_metrics.txt") # get file names
dat7<-readFiles(files7, 7)
files8 <- dir(data_path, pattern = "*8_metrics.txt") # get file names
dat8<-readFiles(files8, 8)
files9 <- dir(data_path, pattern = "*9_metrics.txt") # get file names
dat9<-readFiles(files9, 9)
files10 <- dir(data_path, pattern = "*10_metrics.txt") # get file names
dat10<-readFiles(files10, 10)

####################################################################################
#reload the results if needed
setwd("/Users/zhangxiaopu/Desktop/Isoform/02.modelTraining/02.metrics/")
load("metrics0.RData")
datR<-data2
load("metrics1.RData")
dat1<-data2
load("metrics2.RData")
dat2<-data2
load("metrics3.RData")
dat3<-data2
load("metrics4.RData")
dat4<-data2
load("metrics5.RData")
dat5<-data2
load("metrics6.RData")
dat6<-data2
load("metrics7.RData")
dat7<-data2
load("metrics8.RData")
dat8<-data2
load("metrics9.RData")
dat9<-data2
load("metrics10.RData")
dat10<-data2

data3<-rbind(datR,dat1,dat2,dat3,dat4,dat5,dat6,dat7,dat8,dat9,dat10)
breaks<-seq(0,1,0.1)
data3$pvalue_bin<-cut(data3$AccuracyPValue, breaks=breaks, include.lowest=TRUE, right=FALSE)

pvalue_counts<-data3 %>% group_by(pvalue_bin,model, perm) %>% tally()
pvalue_counts<-pvalue_counts %>% mutate_if(is.numeric, funs(replace_na(., 0)))
pvalue_counts<-pvalue_counts %>% spread(perm, n)


pvalue_counts$ratio1<-(pvalue_counts$metrics0+1)/(pvalue_counts$`metrics1`+1)
pvalue_counts$ratio2<-(pvalue_counts$metrics0+1)/(pvalue_counts$`metrics2`+1)
pvalue_counts$ratio3<-(pvalue_counts$metrics0+1)/(pvalue_counts$`metrics3`+1)
pvalue_counts$ratio4<-(pvalue_counts$metrics0+1)/(pvalue_counts$`metrics4`+1)
pvalue_counts$ratio5<-(pvalue_counts$metrics0+1)/(pvalue_counts$`metrics5`+1)
pvalue_counts$ratio6<-(pvalue_counts$metrics0+1)/(pvalue_counts$`metrics6`+1)
pvalue_counts$ratio7<-(pvalue_counts$metrics0+1)/(pvalue_counts$`metrics7`+1)
pvalue_counts$ratio8<-(pvalue_counts$metrics0+1)/(pvalue_counts$`metrics8`+1)
pvalue_counts$ratio9<-(pvalue_counts$metrics0+1)/(pvalue_counts$`metrics9`+1)
pvalue_counts$ratio10<-(pvalue_counts$metrics0+1)/(pvalue_counts$`metrics10`+1)

pvalue <- pvalue_counts %>% select(pvalue_bin, model, ratio1, ratio2, ratio3, ratio4,
                                   ratio5, ratio6, ratio7, ratio8, ratio9, ratio10)
pvalue <- melt(pvalue, id.vars = c("pvalue_bin", "model"))
colnames(pvalue) <- c("pvalue bin", "model", "variable", "Num Real / Num Perm")
pvalue$model[pvalue$model=="1"] <- "logistic regression"
pvalue$model[pvalue$model=="2"] <- "elastic net"
pvalue$model[pvalue$model=="3"] <- "random forest"
pvalue$model[pvalue$model=="4"] <- "xgboost"

p2 <- ggboxplot(pvalue, x="pvalue bin", y="Num Real / Num Perm", color="model") +
  scale_y_continuous(trans='log2')


win <- read.csv("~/Desktop/Isoform/01.winComparison/Gene.Window", header=T, sep = "\t")
win <- melt(win, id.vars = c("window.region", "GeneID"))
win$`log(win num)` <- log(win$value)
colnames(win) <- c("region","ID", "window size", "value", "log(number of windows)")

p1 <- ggboxplot(win, x="window size", y="log(number of windows)", color="window size")

pdf("/Users/zhangxiaopu/Desktop/Isoform/02.modelTraining/Fig2.pdf")
ggarrange(p1, p2, ncol=1, heights=c(1,1.5), labels=c("A","B"),
          font.label = list(size = 14, color = "black", face = "bold", family = NULL))
dev.off()
#this bit gets corrected p value given permutation results
#when have more permutations would compare to results from all not just dat1
#count how many permutation p values were less than each real one for each model type
#correct P value is FDR
datALL <- rbind(dat1,dat2,dat3,dat4,dat5,dat6,dat7,dat8,dat9,dat10)
datR$perm_count<-NA
datR$corrected_p<-NA
for(i in as.numeric(unique(datR$model)))
{
  datR$perm_count[which(datR$model == i)]<- sapply(datR$AccuracyPValue[which(datR$model == i)], function(x) sum(x >= datALL$AccuracyPValue[which(datALL$model == i)]))
  datR$corrected_p[which(datR$model == i)]<-datR$perm_count[which(datR$model == i)]/length(datALL$AccuracyPValue[which(datALL$model == i)])
}
# write.table(datR,"~/correctP.permutation",row.names = F)

################################################################################################
##expression level
                                                   
library(tidyverse)
library(reshape2)
library(ggpubr)

setwd("~/Dropbox/IsoformUsage/SigGenes/")

load("metrics0.RData")
datR<-data2
load("metrics1.RData")
dat1<-data2
load("metrics2.RData")
dat2<-data2
load("metrics3.RData")
dat3<-data2
load("metrics4.RData")
dat4<-data2
load("metrics5.RData")
dat5<-data2
load("metrics6.RData")
dat6<-data2
load("metrics7.RData")
dat7<-data2
load("metrics8.RData")
dat8<-data2
load("metrics9.RData")
dat9<-data2
load("metrics10.RData")
dat10<-data2

data3<-rbind(datR,dat1,dat2,dat3,dat4,dat5,dat6,dat7,dat8,dat9,dat10)
breaks<-seq(0,1,0.1)
data3$pvalue_bin<-cut(data3$AccuracyPValue, breaks=breaks, include.lowest=TRUE, right=FALSE)
data3$pvalue_bin<-str_replace_all(data3$pvalue_bin, "[\\[\\)\\]]", "")
data3$pvalue_bin<-str_replace_all(data3$pvalue_bin, "\\,", "-")


expression<-read_tsv("allRPKM")
data4 <- merge(data3, expression, by.x="gene")
#get the approximate values that splits it into thirds
quantile(data4$RPKM_EUR, c(.3333, .666))

data4$expression_level[data4$RPKM_EUR <= 0.005] <- "low"
data4$expression_level[data4$RPKM_EUR >= 1] <- "high"
data4$expression_level[data4$RPKM_EUR > 0.005 & data4$RPKM_EUR < 1] <- "medium"


pvalue_counts<- data4 %>% filter(model %in% c(2,4)) %>% 
  group_by(pvalue_bin,model, perm, expression_level) %>% tally()
pvalue_counts<-pvalue_counts %>% spread(perm, n)
pvalue_counts<-pvalue_counts %>% mutate_if(is.numeric, funs(replace_na(., 0)))

pvalue_counts$ratio1<-(pvalue_counts$metrics0+1)/(pvalue_counts$`metrics1`+1)
pvalue_counts$ratio2<-(pvalue_counts$metrics0+1)/(pvalue_counts$`metrics2`+1)
pvalue_counts$ratio3<-(pvalue_counts$metrics0+1)/(pvalue_counts$`metrics3`+1)
pvalue_counts$ratio4<-(pvalue_counts$metrics0+1)/(pvalue_counts$`metrics4`+1)
pvalue_counts$ratio5<-(pvalue_counts$metrics0+1)/(pvalue_counts$`metrics5`+1)
pvalue_counts$ratio6<-(pvalue_counts$metrics0+1)/(pvalue_counts$`metrics6`+1)
pvalue_counts$ratio7<-(pvalue_counts$metrics0+1)/(pvalue_counts$`metrics7`+1)
pvalue_counts$ratio8<-(pvalue_counts$metrics0+1)/(pvalue_counts$`metrics8`+1)
pvalue_counts$ratio9<-(pvalue_counts$metrics0+1)/(pvalue_counts$`metrics9`+1)
pvalue_counts$ratio10<-(pvalue_counts$metrics0+1)/(pvalue_counts$`metrics10`+1)

pvalue <- pvalue_counts %>% select(pvalue_bin, model, expression_level, ratio1, ratio2, ratio3, ratio4,
                                   ratio5, ratio6, ratio7, ratio8, ratio9, ratio10)
pvalue <- melt(pvalue, id.vars = c("pvalue_bin", "model","expression_level"))
colnames(pvalue) <- c("pvalue bin", "model","expression", "variable", "Num Real / Num Perm")
pvalue$model[pvalue$model=="2"] <- "elastic net"
pvalue$model[pvalue$model=="4"] <- "xgboost"

pvalueElastic <- pvalue %>% filter(model == "elastic net")
pvalueXgb <- pvalue %>% filter(model=="xgboost")

pvalueElastic$expression <- factor(pvalueElastic$expression, levels = c("high", "medium", "low"))
pvalueXgb$expression <- factor(pvalueXgb$expression, levels = c("high", "medium", "low"))

p4 <- ggboxplot(pvalueElastic, x="pvalue bin", y="Num Real / Num Perm", color="expression") +
  scale_y_continuous(trans='log10', limits=c(0.15, 1000), breaks=c(1,10,100,1000))+ 
  geom_hline(yintercept = 1, lty=2, colour="grey")+ xlab("P value bin") +
  annotation_logticks(sides="l")
p5 <- ggboxplot(pvalueXgb, x="pvalue bin", y="Num Real / Num Perm", color="expression") +
  scale_y_continuous(trans='log10', limits=c(0.15, 1000), breaks=c(1,10,100,1000))+ 
  geom_hline(yintercept = 1, lty=2, colour="grey")+ xlab("P value bin") +
  annotation_logticks(sides="l")

pdf("./FalsePositiveGenes.pdf", width = 15, height = 6)
ggarrange(p4,p5,ncol=2,labels=c("A","B"))
dev.off()
