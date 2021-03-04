library(tidyverse)
library(ggpubr)
library(RColorBrewer)

path_elastic<-"Isoform/03.ModelTraining/02.varImp"
files<-list.files(path=path_elastic,pattern = "xgbTree.txt$", recursive = TRUE) # recursive means listing files sub-dir

sigGenes<-read_delim("Isoform/04.GEUVADIS/xgbTree.sigGene", delim=" ", col_names = FALSE)
sigGenes$name <- paste(sigGenes$X1,"_varImp_elasticnet.txt",sep="")

filename <- intersect(sigGenes$name, files)
all<-list()
for(file in filename)
{
  dat<-read_delim(paste(path_elastic,file, sep=""),col_names = c("Overall", "ID", "Val"), skip=1, delim = " ")
  dat$Overall <- file
  dat$ID <- gsub("exon","exon_",dat$ID)
  dat$ID <- gsub("intron","intron_",dat$ID)
  dat$ID <- gsub("flank","flank_",dat$ID)
  dat <- dat %>% separate(ID, "_", into=c("Region", "RegionNum", "Win")) %>% 
    separate(Overall, "_", into = c("Gene", NA, NA))
  all[[file]]<-dat 
}

varImp<-bind_rows(all)

# exclude genes with maxImp==0
genesWithImp<-varImp %>% group_by(Gene) %>% summarise(maxImp=max(Val)) %>% filter(maxImp > 0)
varImp<-varImp %>% filter(Gene %in% genesWithImp$Gene)
varImp$RegionNum<-as.numeric(varImp$RegionNum)

#convert relative importance to ntiles
varImp<-varImp %>% group_by(Gene, Region) %>% 
  mutate(maxRegionNum = max(RegionNum), fromEnd=(RegionNum-maxRegionNum))

#first randomise rows so variables with importance of 0 arent skewed by window
#only rank windows with a non-zero importance?
varImp<-varImp[sample(nrow(varImp)),]
varImp<-varImp %>% #filter(Val > 0) %>%
  group_by(Gene) %>% mutate(decile_rank = ntile(Val,10))
varImp$Region[which(varImp$Region == "flank" & varImp$RegionNum == 1)]<-"5'Flank"
varImp$Region[which(varImp$Region == "flank" & varImp$RegionNum == 2)]<-"3'Flank"

#this to restrict to multi-exonic genes
exonCounts<-varImp %>% filter(Region=="exon") %>% group_by(Gene) %>% summarise(numExon=max(maxRegionNum)) %>% filter(numExon > 1)
varImp<-varImp %>% filter(Gene %in% exonCounts$Gene[which(exonCounts$numExon >2)])

colours<-brewer.pal(n = 8, name = "Set1")

######5 PRIME PLOT
medians<-varImp %>% ungroup() %>% group_by(Region, RegionNum, Win) %>% summarise(mean=mean(decile_rank), n=n())
medians <- medians %>% arrange(RegionNum,Region, Win)

medians$Label5<-paste(medians$Region, medians$RegionNum," ",medians$Win, sep="")
medians<-medians %>% filter(RegionNum <9, RegionNum > 0, Region!="3'Flank")
medians$Label5 <- factor(medians$Label5, levels = medians$Label5)
#medians %>% filter(Region=="5'Flank" | Region =="3'Flank")
# p1<-ggplot(medians, aes(Label5, mean, colour=Region, shape=Win))+geom_point()+
#   theme(axis.text.x = element_text(angle = 90))+theme_pubr()
# 

df2 <- medians %>% mutate(Xn=as.numeric(Label5))
p3<-ggplot(df2) +
  geom_rect(aes(xmin=Xn-.5, xmax=Xn+.5, ymin=4.5, ymax=7, fill = Region), alpha=0.1, stat="identity") +
  geom_point(aes(x = Xn, y = mean,colour=Region, shape=Win), size = 2.5) + 
  scale_x_continuous(breaks=df2$Xn, labels=df2$Label5,expand = c(0, 0))+
  theme_pubr()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_y_continuous(expand = c(0, 0))+
  scale_shape_manual(name="Region window",
                     values=c(15,16,9),
                     labels=c("1st", "2nd", "3rd"))+
  scale_fill_manual(values = alpha(c(colours[1], colours[2], "white"), 0.2),name="Gene region",
                    breaks=c("5'Flank", "exon", "intron"),
                    labels=c("Flank", "Exon", "Intron"))+
  scale_color_manual(values=c(colours[1], colours[2], colours[3]),name="Gene region",
                     breaks=c("5'Flank", "exon", "intron"),
                     labels=c("Flank", "Exon", "Intron"))+
  geom_hline(yintercept=5.5, linetype="dashed")+ylab("Relative importance (mean decile rank)")+
  ggtitle("5' end of genes")

p5<-ggplot(medians, aes(Label5, n)) + geom_histogram(stat="identity", fill=colours[4], colour="black")+theme_pubr()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+ylab("Number\nof\ngenes")

######3 PRIME PLOTS
medians<-varImp %>% ungroup() %>% group_by(Region, fromEnd, Win) %>% summarise(mean=mean(decile_rank), n=n())
medians <- medians %>% arrange(fromEnd,desc(Region), Win)
medians$Label3<-paste(medians$Region, medians$fromEnd, medians$Win, sep=" ")
medians<-medians %>% filter(fromEnd >-9, fromEnd < 1, Region!="5'Flank")
medians$Label3 <- factor(medians$Label3, levels = medians$Label3)
medians %>% filter(Region=="5'Flank" | Region =="3'Flank")
p2<-ggplot(medians, aes(Label3, mean, colour=Region, shape=Win))+geom_point()+
  theme(axis.text.x = element_text(angle = 90))+ylim(c(7,9.5))+theme_pubr()

df2 <- medians %>% mutate(Xn=as.numeric(Label3))
p4<-ggplot(df2) +
  geom_rect(aes(xmin=Xn-.5, xmax=Xn+.5, ymin=4.5, ymax=7, fill = Region), alpha=0.1, stat="identity") +
  geom_point(aes(x = Xn, y = mean,colour=Region, shape=Win), size = 2.5) + 
  scale_x_continuous(breaks=df2$Xn, labels=df2$Label3,expand = c(0, 0))+
  theme_pubr()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  scale_y_continuous(expand = c(0, 0))+
  scale_shape_manual(name="Region window",
                     values=c(15,16,9),
                     labels=c("1st", "2nd", "3rd", "4th", "5th"))+
  scale_fill_manual(values = alpha(c(colours[1], colours[2], "white"), 0.2))+
  scale_color_manual(values=c(colours[1], colours[2], colours[3]))+
  geom_hline(yintercept=5.5, linetype="dashed")+ggtitle("3' end of genes")

p6<-ggplot(medians, aes(Label3, n)) + geom_histogram(stat="identity", fill=colours[4], colour="black")+theme_pubr()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

##Main plot
#ggarrange(p1,p2, common.legend = TRUE)
pdf("Isoform/03.ModelTraining/varImp_xgb.pdf", width=11, height=6)
ggarrange(p3,p4,p5,p6, common.legend = TRUE, nrow=2, ncol=2, align="hv", heights=c(5,1))
dev.off()


