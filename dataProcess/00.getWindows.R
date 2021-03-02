## generate mapping models
## outputfile should be like:
#chr1 start end gene_exon1 gene
#chr1 start end gene_exon2 gene
#chr1 start end gene_exon3 gene

# step 1: get intron and flank
library(reshape2)
args = commandArgs(trailingOnly=TRUE)
win <- as.numeric(args[1])
flank <- as.numeric(args[2])
outputname <- paste(win,flank,sep = ".")


getWindow <- function(df, win){
  length <- (df$end - df$start) %/% win
  colnames(df) <- c("gene", "win.start1", paste("win.end",win,sep = ""), "chr", "region")
  name <- colnames(df)
    
  for (i in c(1:(win-1))){
    df <- cbind(df, 
                df$win.start1 + length*i, 
                df$win.start1 + length*i +1) ## winstart, end, end+1
    name <- c(name, paste("win.end",i,sep=""), paste("win.start",i+1,sep=""))
  }
  
  ### when a certain region length =0
  df[which(length==0),] <- cbind(df[which(length==0),][,1:5],
                                 data.frame(rep(df[which(length==0),][,c(3,2)],4)))
  colnames(df) <- name
  
  outcome <- data.frame()
  for (i in c(1:win)){
    dftemp <- df[,c("gene", paste("win.start",i,sep=""), paste("win.end",i,sep=""), "region", "chr")]
    colnames(dftemp) <- c("gene", "start", "end", "region", "chr")
    dftemp$region <- paste(dftemp$region, i, sep="_win")
    outcome <- rbind(outcome, dftemp)
  }
  outcome <- outcome[order(outcome$region),]
  outcome <- outcome[,c("chr","start","end","region","gene")]
  return(outcome) 
}

getIntron <- function(df, win){
  df$region <- "exon"
  label <- c()
  for(i in c(1:length(df$region))) {label <- c(label,paste("exon",i,sep=""))}
  
  df.flank <- data.frame(df$gene[1],
                         c(df$start[1]-(flank+1), df$end[nrow(df)]+1 ), 
                         c(df$start[1]-1, df$end[nrow(df)]+(flank+1) ),
                         df$chr[1])
  df.flank$region <- "flank"
  for(i in c(1:length(df.flank$region))) {label <- c(label,paste("flank",i,sep=""))}
  colnames(df.flank) <- c("gene","start","end","chr","region")
  res <- rbind(df, df.flank)
  
  if ( nrow(df) > 1){
    df.intron <- data.frame(df$gene[1],
                            df$end[1:nrow(df)-1]+1, df$start[2:nrow(df)]-1,
                            df$chr[1])
    df.intron$region <- "intron"
    for(i in c(1:length(df.intron$region))) {label <- c(label,paste("intron",i,sep=""))}
    colnames(df.intron) <- c("gene","start","end","chr","region")
    df.intron[which(df.intron$start>df.intron$end),]$start <-df.intron[which(df.intron$start>df.intron$end),]$end
    res <- rbind(res, df.intron)
  }
  
  res$region <- label
  window.df <- getWindow(res, win)
  return(window.df)
}


data <- read.csv("~/Desktop/bedtoolsMerge.chr", header=T, sep="\t", stringsAsFactors = F)

for (name in unique(data$gene)){
  temp <- data[data$gene==name,]
  res <- getIntron(temp, win)
  write.table(res, paste("Isoform/01.OwnReference/Mapping.win",outputname,sep = "."), 
              sep="\t", col.names=F, row.names=F, quote=F, append=T)
  }

