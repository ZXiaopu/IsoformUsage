library(ellipse)
library(GGally)
library(RColorBrewer)
library(tidyverse)

#data_path <- "~/Desktop/Project/1. Isoform/02.modelTraining/03.PermutationAnalysis/FUMA/result"
data_path <- "C:\\Users\\jgdpr\\Dropbox\\IsoformUsage\\FUMA\\result"
files <- list.files(data_path, pattern ="GS")
all <- list()

## For each method, select 10 gene sets with the lowest P value as candidate categories
for (file in files) {
  dat <- read_delim(paste(data_path, file, sep="/"), col_names = T, delim = "\t") %>% 
    mutate(ID=paste(Category, GeneSet, sep="_"))#%>%
    #filter(Category %in% c("Curated_gene_sets", "Wikipathways", "GWAScatalog"))
    #for (c in dat$Category) {
    #  dat1 <- dat %>% filter(Category %in% c)
    #  dat1 <- dat1[order(dat1$adjP),]
    #  dat1 <- dat1[1:10,]
    #  #all[[paste(file, c, sep="_")]] <- dat1
  id<-str_replace(file, "GS_", "")
  id<-str_replace(id, ".txt", "")
  
      all[[id]] <- dat
  #}
}  

all["elasticAndxgb"]<-NULL

GS <- bind_rows(all, .id="model")

GS2<-GS %>% select(ID, model, p) %>% distinct() %>% spread(model, p)

ggpairs(-log10(GS2[,-1]+1e-300), method = c("complete.obs", "pearson"), col = "red")+xlab(expression(-~log[10]~("Enrichment P")))+ylab(expression(-~log[10]~("Enrichment P")))

my_colors <- brewer.pal(5, "Set1")

#Panel of correlations
panel.corr <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y, use="complete.obs"), digits=2)
  txt <- paste0("Corr: ", r)
  text(0.5, 0.5, txt, cex = 1.2)
}

#Panel of histograms
panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks
  len <- length(breaks)
  y <- h$counts/max(h$counts)
  rect(breaks[-len], 0, breaks[-1], y, col = my_colors[1])
}

#Panel of scatterplots
panel.scat <- function(x, y){
  points(x,y, pch = 19, cex = 1, col = my_colors[2])
}

pdf("C:\\Users\\jgdpr\\Dropbox\\IsoformUsage\\FUMA\\correlogram.pdf")
#Plot
pairs(-log10(GS2[,-1]+1e-300),  
      lower.panel = panel.scat,
      upper.panel = panel.corr,
      diag.panel = panel.hist,
      labels = c("edgeR","elastic net","GEUVADIS",
                 "log. regr.","rand. forest","grad. boost"),
      gap = 0.3, 
      main=expression(-~log[10]~("Enrichment P")))
dev.off()

