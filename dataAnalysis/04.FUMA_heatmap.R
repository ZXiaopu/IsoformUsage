library(tidyverse)

data_path <- "/Users/zhangxiaopu/Desktop/Isoform/02.modelTraining/03.PermutationAnalysis//FUMA"

files <- list.files(data_path, pattern ="GS")
all <- list()

for (file in files) {
  dat <- read_delim(paste(data_path, file, sep="/"), col_names = T, delim = "\t") %>%
    filter(Category %in% c("Curated_gene_sets", "Wikipathways", "GWAScatalog"))
  for (c in dat$Category) {
    dat1 <- dat %>% filter(Category %in% c)
    dat1 <- dat1[order(dat1$adjP),]
    dat1 <- dat1[1:10,]
    all[[paste(file, c, sep="_")]] <- dat1
  }
}  

GS <- bind_rows(all)
GS <- GS %>% select(Category, GeneSet, adjP)

res <- list()
for (file in files) {
  dat <- read_delim(paste(data_path, file, sep="/"), col_names = T, delim = "\t") %>%
    filter(Category %in% c("Curated_gene_sets", "Wikipathways", "GWAScatalog")) %>%
    filter(GeneSet %in% unique(GS$GeneSet))
  dat$file <- file
  dat <- dat %>% separate(file, "_", into=c(NA, "temp"))%>%
    separate(temp, sep=".txt", into=c("model",NA))
  dat <- dat %>% select(Category, GeneSet, adjP, model)
  res[[file]] <- dat
}

output <- bind_rows(res)
output$`-log10(adjusted P)` <- -log10(output$adjP)
output <- output %>% select(-adjP) %>% spread(key = model, value = `-log10(adjusted P)`) %>%
  # mutate_if(is.numeric, funs(replace_na(., 0))) %>% 
  as.data.frame()

rownames(output) <- output$GeneSet
colnames(output) <- c("Category", "GeneSet", "Elastic Net", "GEUVADIS", 
                      "Logistic Regression", "Random Forest", "Xgboost",
                      "Xgboost, elastic net Intersect")
output <- output %>% arrange("Category", "GeneSet", "GEUVADIS", 
                             "Logistic Regression", "Elastic Net", "Random Forest", "Xgboost",
                             "Xgboost, elastic net Intersect")
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(0, 10, 20, 30), c("white", "cornflowerblue", "yellow", "red"))
mat<-output[,c(3:8)]
row.names(mat)<-output$GeneSet
ha = rowAnnotation(Category = output$Category,
                   col = list(Category = c("GWAScatalog" = "red", "Curated_gene_sets" = "green", "Wikipathways" = "blue")))
pdf("~/Desktop/Isoform/02.modelTraining/03.PermutationAnalysis/FUMA/Fig4.Heatmap.pdf", width=11, height=13)
Heatmap(as.matrix(mat), col = col_fun, 
        width = ncol(input)*unit(12, "mm"),
        row_names_max_width = unit(8,"cm"),
        row_names_gp = grid::gpar(fontsize = 8),
        right_annotation = ha,  name="-log10(adjusted P)")
dev.off()


input <- output[,3:8]
df_ha <- data.frame(output$Category)
rownames(df_ha) <- rownames(output)
colnames(df_ha) <- "Category"
ha <- rowAnnotation(df=df_ha,
                        col = list(Category = c("Curated_gene_sets" = "red", 
                                                "GWAScatalog" = "blue",
                                                "Wikipathways" = "green")))

input <- as.matrix(output[,3:8])
pdf("~/Desktop/Isoform/03.Permutation/FUMA/heatmap.pdf", width=10, height=13)
draw(Heatmap(input, name="-log10(adjusted P)",
             width = ncol(input)*unit(12, "mm"),
             row_names_max_width = 
             row_names_gp = grid::gpar(fontsize = 7.3),
             right_annotation = ha))
dev.off()
     