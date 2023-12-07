library(plyr)
library(magrittr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

pauRange <- args[1] #"25_75"
inputfile <- args[2] #paste("SP1_S2S3E_peaks_count_table_",pauRange,".tab", sep="")
wd <- args[3] #"C:\\RSYNC\\Nabeel\\clip_analysis\\SP1_2021"
target <- args[4] #"SP1_pureclipPeaks"

setwd(wd) 
 
count.table <- type.convert(read.delim(file.path(wd,inputfile), header = F, stringsAsFactors = F))
colnames(count.table) <- c("chr", "start", "end", "utrID","length", "strand", "geneID", "geneName", "peakCount", "PAU", "utrType")
head(count.table)



for(pau in unique(count.table$PAU)){
  #pau <- "Down"
  putr_sub <- unique(count.table[count.table$utrType=="putr" & count.table$PAU == pau, ])
  targeted_putr <- putr_sub[putr_sub$peakCount > 0, 1:8]
  nonTargeted_putr <- putr_sub[putr_sub$peakCount == 0, 1:8]
  
  write.table(targeted_putr, paste("siSP1_PAU_",pau,"_putr_",pauRange,"_targeted_by_",target,".bed",sep=""), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(nonTargeted_putr, paste("siSP1_PAU_",pau,"_putr_",pauRange,"_nonTargeted_by_",target,".bed",sep=""), row.names=F, col.names=F, sep="\t", quote=F)
  
  dutr_sub <- unique(count.table[count.table$utrType=="dutr" & count.table$PAU == pau, ])
  print(pau)
  print(dim(unique(putr_sub)))
  print(dim(unique(dutr_sub)))
  length(unique(putr_sub$geneID))
  merged.table <- merge.data.table(putr_sub, dutr_sub, by=c("geneID"), all.x=T)
  merged.table <- unique(merged.table)
  print(dim(merged.table))
  gi <- merged.table$geneID
  gi_dup <- gi[duplicated(gi)]
  merged.table[merged.table$geneID==gi_dup,]
  
  merged.table <- mutate(merged.table, count_diff=peakCount.y-peakCount.x) %>%
                  mutate(length_diff=length.y-length.x)
  merged.table_targeted <- merged.table[merged.table$count_diff > 0, ]
  merged.table_nonTargeted <- merged.table[merged.table$count_diff == 0, ]
  targeted_dutr <- dutr_sub[dutr_sub$geneID %in% merged.table_targeted$geneID, 1:8]
  nonTargeted_dutr <- dutr_sub[dutr_sub$geneID %in% merged.table_nonTargeted$geneID, 1:8]
  
  write.table(targeted_dutr, paste("siSP1_PAU_",pau,"_dutr_",pauRange,"_targeted_by_",target,".bed",sep=""), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(nonTargeted_dutr, paste("siSP1_PAU_",pau,"_dutr_",pauRange,"_nonTargeted_by_",target,".bed",sep=""), row.names=F, col.names=F, sep="\t", quote=F)
}


