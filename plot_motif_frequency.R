
# plot frequency of motif occurrence after excuting the following commands on bc cluster:

# count_motif_frequency_for_bed_bc.sh

library(ggpubr)
library(ggplot2)
library(ggsignif)
library(sciplot)
library(plyr)

args = commandArgs(trailingOnly=TRUE)

wdir <- args[1] ##"C:/RSYNC/Jingwen/RNAseq/qapa/PAU_analysis_batch2/results_output/GO_results_25_75"

setwd(wdir)
motif <- args[2] ##"DGGDRG"
peak <- args[3] ##"SP1_S2S3E"
plotd <- args[4]

for(range in c("0_100", "25_75")){
   count.table <- NULL
   for(utr in c("dutr", "putr")){
      for(group in c("Down", "NoChange", "Up")){
        #group <- "down"
         target_bed_file <- paste("siSP1_PAU_",group,"_",utr,"_",range,"_targeted_by_",peak,".bed", sep="")
         target_bed <- read.table(target_bed_file, header=F)
         colnames(target_bed) <- c("chr", "start", "end", "utr", "len", "strand", "gene", "name")
         target_bed <- mutate(target_bed, seqID=paste(chr, ":", start, "-", end, "(", strand, ")", sep=""))
         
         target_motif_count_file <- paste("siSP1_PAU_",group,"_",utr,"_",range,"_targeted_by_",peak,"_",motif,"_count.tab", sep="")
         target_motif_count_table <- read.table(target_motif_count_file, header=T)
         dim(target_motif_count_table)
         target_sequence_without_motif <- dim(target_bed)[1] - dim(target_motif_count_table)[1]
         target_motif_count_merged_table <- merge(target_motif_count_table, target_bed, by="seqID")
         
         target_motif_count_merged_table <- mutate(target_motif_count_merged_table, freq_per_kb=(frequency/ceiling(len/1000)))
         
         target_motif_count_perkb <- target_motif_count_merged_table$freq_per_kb
        
         
         non_target_bed_file <- paste("siSP1_PAU_",group,"_",utr,"_",range,"_nonTargeted_by_",peak,".bed", sep="")
         non_target_bed <- read.table(non_target_bed_file, header=F)
         colnames(non_target_bed) <- c("chr", "start", "end", "utr", "len", "strand", "gene", "name")
         non_target_bed <- mutate(non_target_bed, seqID=paste(chr, ":", start, "-", end, "(", strand, ")", sep=""))
         
         non_target_motif_count_file <- paste("siSP1_PAU_",group,"_",utr,"_",range,"_nonTargeted_by_",peak,"_",motif,"_count.tab", sep="")
         non_target_motif_count_table <- read.table(non_target_motif_count_file, header=T)
         dim(non_target_motif_count_table)
         non_target_sequence_without_motif <- dim(non_target_bed)[1] - dim(non_target_motif_count_table)[1]
         non_target_motif_count_merged_table <- merge(non_target_motif_count_table, non_target_bed, by="seqID")
         
         non_target_motif_count_merged_table <- mutate(non_target_motif_count_merged_table, freq_per_kb=(frequency/ceiling(len/1000)))
         
         non_target_motif_count_perkb <- non_target_motif_count_merged_table$freq_per_kb
         
         
         target_motif_count_perkb_complete <- c(rep(0, target_sequence_without_motif), target_motif_count_perkb)
         non_target_motif_count_perkb_complete <- c(rep(0, non_target_sequence_without_motif), non_target_motif_count_perkb)
         
         if(0){ ## old way of plot, not used
         fr <- var.test(target_motif_count_perkb_complete, non_target_motif_count_perkb_complete, alternative = "two.sided")
         print(fr$p.value)
         
         VAR <- TRUE
         if(fr$p.value < 0.05){
           VAR <- FALSE
         }
         
         tr <- t.test(target_motif_count_perkb_complete, non_target_motif_count_perkb_complete, alternative = "two.sided", var.equal = VAR)
         print(tr$p.value)
         
         perkb_list <- list()
         perkb_list[["targeted"]] <- target_motif_count_perkb_complete
         perkb_list[["non-targeted"]] <- non_target_motif_count_perkb_complete
         Group <- names(perkb_list)
         Mean <- unlist(lapply(perkb_list, mean))
         SE <- unlist(lapply(perkb_list, se))
         Df <- data.frame(Group, Mean, SE)
         
         #png(paste("PAU_analysis_siSP1_", group, "_list_dutr_", motif, "_barplot.png", sep=""), width=800, height=800)
         
         mp <- ggplot(Df, aes(x = Group, y = Mean)) +
           geom_bar(stat = "identity", fill="blue") +
           geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = .2) +
           geom_signif(comparisons = list(c("targeted", "non-targeted")), annotations=c(tr$p.value), y_position = c(Mean + SE + 0.1), tip_length = 0.1) +
           ggtitle(paste("Frequency of motif", motif, sep=" ")) + labs(y = "Average # of motifs per kb", x = "") +
           theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(colour = "black", size = rel(1.2)),
                 axis.title = element_text(colour = "black", size = rel(1.5)))
         
         ggsave(filename=paste("siSP1_PAU_",group,"_",utr,"_",range,"_targeted_by_",peak,"_",motif, "_barplot.pdf", sep=""), plot=mp,
                width = 5, height = 5, 
                units = "in", # other options are "in", "cm", "mm" 
                dpi = 200
         )
	 }

	targeted <- as.data.frame(target_motif_count_perkb_complete) %>%
			mutate(TARGET="targeted", UTR=utr, PAU=group)
	colnames(targeted) <- c("motif_frequency", "target", "utr", "pau")
	non_targeted <- as.data.frame(non_target_motif_count_perkb_complete)  %>%
                        mutate(TARGET="non-targeted", UTR=utr, PAU=group)
	colnames(non_targeted) <- c("motif_frequency", "target", "utr", "pau")
	count.table <- rbind(count.table, targeted, non_targeted)

      }
   }
   #count.table <- as.data.frame(count.table)
   #colnames(count.table) <- c("motif_frequency", "target", "utr", "pau")
   p <- ggbarplot(count.table, "target", "motif_frequency",
               color = "target", fill="target", palette = "ucscgb", title=paste("Frequency of motif", motif),
               xlab=FALSE, facet.by=c("pau", "utr"), add="mean_se", ylab="Motif frequency (mean+/-se)") %>%
        ggpar(legend = "none", ylim=c(0,10)) +
        stat_compare_means(label = "p.format", label.x.npc="center", label.y=9.0)

  figname <- file.path(plotd, paste(motif, "_motif_frequency_in_UTR_", range, "_targeted_by_", peak, "_bar.pdf", sep=""))
  ggexport(p, filename=figname)

  p <- ggboxplot(count.table, "target", "motif_frequency", ylab="Motif frequency", 
               color = "target", palette = "ucscgb", title=paste("Frequency of motif", motif),
               xlab=FALSE, facet.by=c("pau", "utr")) %>%
        ggpar(legend = "none") +
        stat_compare_means(label = "p.format", label.x.npc="center", label.y.npc=0.9)

  figname <- file.path(plotd, paste(motif, "_motif_frequency_in_UTR_", range, "_targeted_by_", peak, "_box.pdf", sep=""))
  ggexport(p, filename=figname)

   
  p <- ggviolin(count.table, "target", "motif_frequency", ylab="Motif frequency", 
               color = "target", palette = "ucscgb", title=paste("Frequency of motif", motif),
               xlab=FALSE, facet.by=c("pau", "utr"), add="boxplot") %>%
	ggpar(legend = "none") +
	stat_compare_means(label = "p.format", label.x.npc="center", label.y.npc=0.9)

  figname <- file.path(plotd, paste(motif, "_motif_frequency_in_UTR_", range, "_targeted_by_", peak, "_violin.pdf", sep=""))
  ggexport(p, filename=figname)

}

