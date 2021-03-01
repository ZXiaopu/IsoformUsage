### GEUVADIS file downloaded from https://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/files/analysis_results/?ref=E-GEUV-1

library(tidyverse)
setwd("Isoform/04.GEUVADIS")
trans <- read.csv("GD660.TrQuantRPKM.txt",sep="\t")
gene <- read.csv("GD660.GeneQuantRPKM.txt",sep = "\t")

Down462 <- read.csv("E-GEUV-1.sdrf.txt",sep = "\t", stringsAsFactors = F)
Down462$Extract.Name <- gsub(" extract","",Down462$Extract.Name)
data462 <- Down462 %>% select(Characteristics.ancestry.category.,Extract.Name)
colnames(data462) <- c("ancestry","name")

data462$ancestry[which(data462$ancestry %in% c("British","Tuscan","Utah","Finnish"))] <- "EUR"
data462$ancestry[which(data462$ancestry %in% c("Yoruba"))] <- "YRI"
data462 <- unique(data462)

for(tID in as.character(unique(trans$TargetID)))
{
  transRow <- trans[trans$TargetID==tID,]
  geneRow <- gene[gene$Gene_Symbol==transRow$Gene_Symbol,]
  inf <- geneRow[,c(1:4)]
  trRatio <- transRow[,-c(1:4)]/geneRow[,-c(1:4)]
  trRatio <- trRatio %>% mutate_if(is.numeric, funs(replace_na(., 0)))
  trRatio1 <- data.frame(t(trRatio), stringsAsFactors=F)
  trRatio1$name <- rownames(trRatio1)
  res <- merge(trRatio1, data462, by.x="name")
  colnames(res) <- c("name","ratio","ancestry")
  EUR <- as.numeric(res$ratio[res$ancestry=="EUR"])
  YRI <- as.numeric(res$ratio[res$ancestry=="YRI"])
  res <- wilcox.test(EUR, YRI, paired=FALSE)
  output <- t(c(tID, as.character(transRow$Gene_Symbol[transRow$TargetID==tID]),res$p.value))
  write.table(resRatio,"./GEUVADIS.pval",quote = F, sep="\t",
              col.names = F,row.names = F,append = T)
}


pvalue <- read.csv("GEUVADIS.pval",sep="\t",header = F)
colnames(pvalue) <- c("transID", "geneID", "pval")
pvalue$adjustP <- p.adjust(pvalue$pval, method = "BH", n= length(pvalue$pval))
write.table(pvalue,"GEUVADIS.adjustP",col.names = T, row.names = F)
