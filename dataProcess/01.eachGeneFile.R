input <- commandArgs(trailingOnly = FALSE)
file <- as.character(input[1])

setwd("Isoform/01.OwnReference/MapRes")
popfile <- read.csv("populationFile", header=F,sep = "\t")
sampleID <- gsub(".Map.win3.flank1500","",file)
pop <- popfile$V1[popfile$V2==sampleID]
data <- read.csv(file, header = F, sep = "\t")
colnames(data) <- c("chr","start","end","region","geneID","reads","coverage","rate")
data1 <- data[,c("geneID", "region", "reads")]

setwd("Isoform/02.Data/win3.flank1500")  
for (gene in unique(data1$geneID)) {
  df <- data1[data1$geneID==gene,]
  df1 <- data.frame(t(df), stringsAsFactors = F)
  colnames(df1) <- df$region
  df1 <- df1[3,]
  df1 <- df1 %>% mutate_if(is.character, as.numeric)
  df1 <- data.frame(t(apply(df1, 1, function(x) x/sum(x))))
  df1$pop <- pop
  write.table(df1, paste(gene, ".win3.flank1500", sep = ""),
              quote = F, row.names = F, col.names=T, append = T)
}
