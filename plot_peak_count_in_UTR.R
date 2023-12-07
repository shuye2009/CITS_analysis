
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(sciplot)
library(plyr)
library(magrittr)

se <- function(x, na.rm=T)
{
  if(na.rm){
    return(sd(x,na.rm=T)/sqrt(sum(!is.na(x))))
  }else{
    x[is.na(x)] <- 0
    return(sd(x)/sqrt(length(x)))
  }
}

args = commandArgs(trailingOnly=TRUE)

inputfile <- args[1] ##"SP1_S2S3E_peaks_count_table_0_100.tab"
wd <- args[2] ##"C:\\RSYNC\\Nabeel\\clip_analysis\\SP1_2021"
target <- args[3] ##"SP1_S2S3E"
pauRange <- args[4] ##"0_100"
plotd <- args[5]

setwd(wd)

pseudocount <- 0.1 ## as an alternative to ceiling
  
my.dist.table <- type.convert(read.delim(file.path(wd,inputfile), header = F, stringsAsFactors = F))
colnames(my.dist.table) <- c("chr", "start", "end", "utrID","length", "strand", "geneID", "geneName", "peakCount", "PAU", "utrType")
  
count.table <- my.dist.table %>% mutate(normCount = peakCount/ceiling(length/1000)) ## counts per kb, use ceiling function to avoid divide by very small numbers

head(count.table)

p <- ggbarplot(count.table, "PAU", "normCount",
               color = "PAU", fill="PAU", palette = "ucscgb", title="Peak density in 3'UTR",
               xlab=FALSE, facet.by="utrType", add="mean_se", ylab="Normalized count (mean+/-se)") %>%
  ggpar(legend = "none", ylim=c(0,3)) +
  stat_compare_means(label = "p.format", label.y=3) +
  stat_compare_means(label = "p.signif", label.y=1.6, tip.length=0.01, step.increase=0.02,
                     comparisons=list(c("Up","Down"), c("Down", "NoChange"), c("Up", "NoChange")), hide.ns = FALSE)

figname <- file.path(plotd, paste(target, "_peak_density_in_UTR_", pauRange, "_bar.pdf", sep=""))
ggexport(p, filename=figname)

p <- ggboxplot(count.table, "PAU", "normCount",  ylab="Normalized count",
               color = "PAU", palette = "ucscgb", title="Peak density in 3'UTR",
               xlab=FALSE, facet.by="utrType") %>%
  ggpar(legend = "none") +
  stat_compare_means(label = "p.format", label.y=12.5) +
  stat_compare_means(label = "p.signif",
                     comparisons=list(c("Up","Down"), c("Down", "NoChange"), c("Up", "NoChange")), hide.ns = FALSE)

figname <- file.path(plotd, paste(target, "_peak_density_in_UTR_", pauRange, "_box.pdf", sep=""))
ggexport(p, filename=figname)

p <- ggviolin(count.table, "PAU", "normCount",  ylab="Normalized count",
               color = "PAU", palette = "ucscgb", title="Peak density in 3'UTR",
               xlab=FALSE, facet.by="utrType", add="boxplot") %>%
  ggpar(legend = "none") +
  stat_compare_means(label = "p.format", label.y=12.5) +
  stat_compare_means(label = "p.signif",
                     comparisons=list(c("Up","Down"), c("Down", "NoChange"), c("Up", "NoChange")), hide.ns = FALSE)

figname <- file.path(paste(target, "_peak_density_in_UTR_", pauRange, "_violin.pdf", sep=""))
ggexport(p, filename=figname)


## THE END


