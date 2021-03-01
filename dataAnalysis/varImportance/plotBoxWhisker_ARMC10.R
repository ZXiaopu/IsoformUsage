library(tidyverse)
library(RColorBrewer)
library(ggpubr)
setwd("~/Desktop/Project/1. Isoform/02.modelTraining/04.varImp/ARMC10/")
data1 <- read.csv("ENSG00000170632.win3.flank1500",sep = "\t",
                  check.names = FALSE,header = TRUE,stringsAsFactors = FALSE)
colnames(data1) <- gsub("exon","exon_",colnames(data1))
colnames(data1) <- gsub("intron","intron_",colnames(data1))
colnames(data1) <- gsub("flank","flank_",colnames(data1))
data1$pop[which(data1$pop %in% c("British", "Finnish", "Tuscan", "Utah"))]<-"EUR"
data1$pop[which(data1$pop %in% c("Yoruba"))]<-"YRI"

len <- read.csv("ENSG00000170632.len",sep = "\t",header = FALSE, 
                col.names = c("chr", "start", "end", "id", "gene"))
len$length <- (len$end-len$start)+1
len$id <- gsub("exon", "exon_", len$id)
len$id <- gsub("intron", "intron_", len$id)
len$id <- gsub("flank", "flank_", len$id)
len <- len %>% select(gene, id, length)

#drop individuals with less than ten reads for now
data1<-data1[which(rowSums(data1[,-c(dim(data1)[2])]) > 10),]
data1[,-c(dim(data1)[2])] <- data1[,-c(dim(data1)[2])]/rowSums(data1[,-c(dim(data1)[2])])

data1 <- gather(data1, "id", "value", -pop) %>% as_tibble()
data1<-left_join(data1,len)
data1$normVal<-data1$value/data1$length

data1 <- data1 %>% separate(id, sep = "_", into=c("Region", "RegionNum", "Window"), convert=TRUE)

data1$label<-paste(data1$Region,data1$RegionNum," ", data1$Window, sep="")

result <- data1 %>% group_by(gene, Region, Window, pop) %>% summarise(median=median(normVal)) %>% 
  spread(pop, median)

data1<-data1 %>% arrange(RegionNum, Region, Window)
nonFlank<-data1 %>% filter(Region!= "flank")
lev<-c("flank1 win1", "flank1 win2", "flank1 win3", unique(nonFlank$label), "flank2 win1", "flank2 win2", "flank2 win3", "flank2 win4", "flank2 win5")
data1$label <- factor(data1$label, levels = lev)

data1 <- data1 %>% mutate(Xn=as.numeric(label))
colnames(data1) <- c("population", colnames(data1)[2:14])

colours<-brewer.pal(n = 9, name = "Set1")
data2<-data1 %>% select(Region, RegionNum, Window, Xn) %>% unique()


p1<-ggplot() + geom_boxplot(data=data1, aes(as.factor(Xn), log10(normVal), colour=population)) +
  geom_rect(data=data2, aes(xmin=Xn-.5, xmax=Xn+.5, ymin=-7, ymax=-3, fill = Region), alpha=0.1, stat="identity") + 
  coord_flip()+theme_pubr()+
  scale_fill_manual(values = alpha(c(colours[2], colours[1], "white"), 0.2),name="Gene region")+
  scale_colour_manual(values = c(colours[3], colours[4]),name="Population")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank(), 
        axis.text=element_text(size=15),
        legend.text = element_text(size=15),
        axis.title.x = element_text(size = 16),
        legend.key = element_rect(fill = "white", colour = "black"))+ 
  annotate("text", y = -7.2, x = 1, label = "5' flank")+
  annotate("text", y = -7.2, x = 63, label = "3' flank")+
  ylab("Log10(Read proportion/window length)")

data3<-data1 %>% filter(Region=="exon", RegionNum==10)
p2<-ggplot() + geom_boxplot(data=data3, aes(label, log10(normVal), colour=population))+coord_flip()+theme_pubr()+
  scale_colour_manual(values = c(colours[3], colours[4]),name="Population")+
  theme(axis.title.y=element_blank(),
        axis.text=element_text(size=30),
        legend.text = element_text(size=30),
        legend.title = element_text(size=30),
        axis.title.x = element_text(size = 30),
        legend.key = element_rect(fill = "white", colour = "black"))+
  ylab("Log10(Read proportion/window length)")
data3$`log10(normVal)` <- log10(data3$normVal)
p2_1 <- ggboxplot(data3,x="label",y="log10(normVal)",color="population")+
  stat_compare_means(aes(group=population),method = "t.test",label="p.signif",size=10)+
  coord_flip()+theme_pubr()+
  scale_colour_manual(values = c(colours[3], colours[4]),name="Population")+
  theme(axis.title.y=element_blank(),
        axis.text=element_text(size=30),
        legend.text = element_text(size=30),
        legend.title = element_text(size=30),
        axis.title.x = element_text(size = 30),
        legend.key = element_rect(fill = "white", colour = "black"))+
  ylab("Log10(Read proportion/window length)")
data4<-read_delim("ENSG00000170632_varImp_elasticnet.txt",
                  col_names = c("Overall", "ID", "Val"), skip=1, delim = " ")
data4$ID <- gsub("exon","exon_",data4$ID)
data4$ID <- gsub("intron","intron_",data4$ID)
data4$ID <- gsub("flank","flank_",data4$ID)

data4 <- data4 %>% separate(ID, sep = "_", into=c("Region", "RegionNum", "Win")) %>%
  filter(Val > 0.1) %>% arrange(Val)
data4$label<-paste(data4$Region,data4$RegionNum, " ", data4$Win, sep="")
data4$label <- factor(data4$label, levels = data4$label)
data4$highlight<-0
data4$highlight[which(data4$Region=="exon" & data4$RegionNum==10 )]<-1
p3<-ggplot(data4, aes(label, Val))+geom_point(aes(colour=as.factor(highlight)), size=4)+ 
  geom_segment( aes(x=label, xend=label, y=0, yend=Val))+coord_flip()+theme_pubr()+ylab("Variable importance (elastic net)")+
  theme(axis.title.y=element_blank(), 
        axis.text=element_text(size=30),
        legend.text = element_text(size=30),
        axis.title.x = element_text(size = 30),
        legend.position = "none")+scale_colour_manual(values = c(colours[9], colours[5]))


data5<-read_delim("ENSG00000170632_varImp_xgbTree.txt",
                  col_names = c("Overall", "ID", "Val"), skip=1, delim = " ")
data5$ID <- gsub("exon","exon_",data5$ID)
data5$ID <- gsub("intron","intron_",data5$ID)
data5$ID <- gsub("flank","flank_",data5$ID)

data5 <- data5 %>% separate(ID, sep = "_", into=c("Region", "RegionNum", "Win")) %>%
  filter(Val > 5) %>% arrange(Val)
data5$label<-paste(data5$Region,data5$RegionNum, " ", data5$Win, sep="")
data5$label <- factor(data5$label, levels = data5$label)
data5$highlight<-0
data5$highlight[which(data5$Region=="exon" & data5$RegionNum==10 )]<-1
p4<-ggplot(data5, aes(label, Val))+geom_point(aes(colour=as.factor(highlight)), size=4)+ 
  geom_segment( aes(x=label, xend=label, y=0, yend=Val))+coord_flip()+theme_pubr()+ylab("Variable importance (xgboost)")+
  theme(axis.title.y=element_blank(), 
        axis.text=element_text(size=30),
        legend.text = element_text(size=30),
        axis.title.x = element_text(size = 30),
        legend.position = "none")+scale_colour_manual(values = c(colours[9], colours[5]))

pBCD <- ggarrange(p2,"",
                  ggarrange(p3,p4,nrow=1,
                            labels=c("C","D"),font.label = list(size = 30)),
                  nrow=1,widths= c(8,1,15), labels=c("B"),font.label = list(size = 30))
pdf("./presentation.pdf", height=20, width = 25)
pBCD1 <- ggarrange(ggarrange(p3,p4,ncol=1,heights = c(0.8,1), align=c("h"),vjust=1,
                             labels=c("A","B"),font.label = list(size = 30)),
                   p2_1,
                   ncol=1,heights= c(2,1), labels=c("","C"), align=c("h"),vjust=1,font.label = list(size = 30))

ggarrange(pBCD1, "", ncol=2, widths = c(0.6,1), labels=c("","D"), font.label = list(size = 30))
dev.off()

pdf("~/Desktop/ARMC10_boxplot.pdf", height=15,width=22)
ggarrange(p1,
          ggarrange(pBCD, "", nrow=2, heights = c(0.6,1), labels=c("","E"), font.label = list(size = 20)),
          widths=c(0.5,1), labels=c("A"), font.label = list(size = 20))

dev.off()

