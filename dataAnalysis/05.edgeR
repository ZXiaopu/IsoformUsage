##########################################################################################  
####################################### run edgeR ########################################
library(tidyverse)
library(edgeR)
setwd("Isoform/05.edgeR")
sample <- read.csv("./sampleTable", row.names = 1,head = TRUE, sep = "\t")
Group <- factor(sample$condition)
design <- model.matrix(~Group)

files <- list.files("./output/",pattern=".txt$") ## files processed by HTseq
myfiles <- lapply(files, function(x) read_delim(paste("./output/",x,sep=""), col_names = c("ID",gsub(".output.txt","",x)),delim="\t"))
myfiles1 <- bind_cols(myfiles)
Counts <- myfiles1[,colnames(myfiles1) %in% rownames(sample)]

exon <- gsub(":",".", myfiles[[1]]$ID)
genename <- substr(exon,1,15)
genes <- data.frame(GeneID=genename, Gene.Exon=exon)

dge <- DGEList(counts=Counts, genes=genes)
dge1 <- dge[filterByExpr(dge, design),]
dgefit <- glmFit(dge1, design, dispersion=0.05)

ds <- diffSpliceDGE(dgefit, geneid="GeneID")
save(ds, "edgeR.RData")


exonLevel <- cbind(ds$genes, ds$exon.p.value)
geneLevel <- cbind(ds$gene.genes, ds$gene.p.value, ds$gene.Simes.p.value)

geneLevel$adjust.Ftest <- p.adjust(geneLevel$`ds$gene.p.value`, method="BH")
geneLevel$adjust.Simes <- p.adjust(geneLevel$`ds$gene.Simes.p.value`, method="BH")

geneLevel <- read_delim("./edgeR.genePvalue", delim = "\t", col_names = T)

write.table(exonLevel, "edgeR.exonPvalue",quote = F, col.names = F, row.names = F)
write.table(geneLevel, "edgeR.genePvalue",quote = F, col.names = F, row.names = F)

sigGenes <- geneLevel %>% filter(adjust.Ftest <=0.05 & adjust.Simes <=0.05)
write.table(data.frame(sigGenes$GeneID),"../02.metrics/edgeR.sigGene", 
            quote = F, col.names = F, row.names = F)
