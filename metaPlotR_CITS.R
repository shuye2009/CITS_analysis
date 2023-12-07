
library("scales")
library("ggplot2")

args = commandArgs(trailingOnly=TRUE)

fact <- args[1]
wd <- args[2] ## "/home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/SP1_analysis/mergedBW"

#setwd("c:/RSYNC/Nabeel/ZBTB48_M6A/results/metaPlotR_plot")
setwd(wd)

inputFile1 <- paste(fact, "_over_Input_DESeq2_insignificant_crossLink_site.bed.dist.measures.txt", sep="")
inputFile2 <- paste(fact, "_over_Input_DESeq2_significant_crossLink_site.bed.dist.measures.txt", sep="")


m6a.dist <- read.delim (inputFile1, header = T)

# Determine longest length transcript for each gene
trx_len <- m6a.dist$utr5_size + m6a.dist$cds_size + m6a.dist$utr3_size # Determine transcript length
temp <- data.frame(m6a.dist$gene_name, m6a.dist$refseqID, trx_len)
colnames(temp) <- c("gene_name", "gid", "trx_len") 
temp.df <- temp[order(temp$gene_name, -temp$trx_len),]
temp.df <- temp.df[!duplicated(temp.df$gene_name),]
temp.df <- temp.df[order(temp.df$gene_name),]

#length(temp.df$gene_name)
#length(temp.df$gene_name[duplicated(temp.df$gene_name)])
#length(temp.df$gene_name[!duplicated(temp.df$gene_name)])

# limit m6a data to one transcript per gene (longest)
m6a.dist <- m6a.dist[m6a.dist$refseqID %in% temp.df$gid,]

# View size of our dataset (rows, columns)
print("nonsignificant peaks table")
dim(m6a.dist)

#qplot(m6a.dist$rel_location, geom="histogram") + geom_vline(xintercept = 1:2, col = "grey") + theme_bw()

#summary(data.frame(m6a.dist$utr5_size, m6a.dist$cds_size, m6a.dist$utr3_size))
utr5.SF <- median(m6a.dist$utr5_size, na.rm = T)/median(m6a.dist$cds_size, na.rm = T)
utr3.SF <- median(m6a.dist$utr3_size, na.rm = T)/median(m6a.dist$cds_size, na.rm = T)

utr5.m6a.dist <- m6a.dist[m6a.dist$rel_location < 1, ]
cds.m6a.dist <- m6a.dist [m6a.dist$rel_location < 2 & m6a.dist$rel_location >= 1, ]
utr3.m6a.dist <- m6a.dist[m6a.dist$rel_location >= 2, ]


utr5.m6a.dist$rel_location <- rescale(utr5.m6a.dist$rel_location, to = c(1-utr5.SF, 1), from = c(0,1))
utr3.m6a.dist$rel_location <- rescale(utr3.m6a.dist$rel_location, to = c(2, 2+utr3.SF), from = c(2,3))

m6a.metagene.coord <- c(utr5.m6a.dist$rel_location, cds.m6a.dist$rel_location, utr3.m6a.dist$rel_location)
#p <- qplot(m6a.metagene.coord, geom="histogram") + geom_vline(xintercept = 1:2, col = "grey") + theme_bw()


#qplot(m6a.metagene.coord, geom="freqpoly") + geom_vline(xintercept = 1:2, col = "red") + theme_bw()

#qplot(m6a.metagene.coord, geom="density") + geom_vline(xintercept = 1:2, col = "red") + theme_bw()

#qplot(m6a.dist$utr3_st, geom="histogram") + xlim (-500,500) + geom_vline(xintercept = 1:2, col = "grey") + theme_bw()


## plot multiple sets of data
#pseudoU.dist <- read.delim("NUDT21_over_Input_DESeq2_significant_crossLink_site.bed.dist.measures.txt", header = T)
pseudoU.dist <- read.delim(inputFile2, header = T)

# Determine longest length transcript for each gene
trx_len <- pseudoU.dist$utr5_size + pseudoU.dist$cds_size + pseudoU.dist$utr3_size # Determine transcript length
temp <- data.frame(pseudoU.dist$gene_name, pseudoU.dist$refseqID, trx_len)
colnames(temp) <- c("gene_name", "gid", "trx_len") 
temp.df <- temp[order(temp$gene_name, -temp$trx_len),]
temp.df <- temp.df[!duplicated(temp.df$gene_name),]
temp.df <- temp.df[order(temp.df$gene_name),]

#length(temp.df$gene_name)
#length(temp.df$gene_name[duplicated(temp.df$gene_name)])
#length(temp.df$gene_name[!duplicated(temp.df$gene_name)])

# limit pseudoU data to one transcript per gene (longest)
pseudoU.dist <- pseudoU.dist[pseudoU.dist$refseqID %in% temp.df$gid,]

print("significant peak table")
dim(pseudoU.dist)

utr5.pseudoU.dist <- pseudoU.dist[pseudoU.dist$rel_location < 1, ]
cds.pseudoU.dist <- pseudoU.dist [pseudoU.dist$rel_location < 2 & pseudoU.dist$rel_location >= 1, ]
utr3.pseudoU.dist <- pseudoU.dist[pseudoU.dist$rel_location >= 2, ]

utr5.SF <- median(pseudoU.dist$utr5_size, na.rm = T)/median(pseudoU.dist$cds_size, na.rm = T)
utr3.SF <- median(pseudoU.dist$utr3_size, na.rm = T)/median(pseudoU.dist$cds_size, na.rm = T)

utr5.pseudoU.dist$rel_location <- rescale(utr5.pseudoU.dist$rel_location, to = c(1-utr5.SF, 1), from = c(0,1))
utr3.pseudoU.dist$rel_location <- rescale(utr3.pseudoU.dist$rel_location, to = c(2, 2+utr3.SF), from = c(2,3))


pseudoU.metagene.coord <- c(utr5.pseudoU.dist$rel_location, cds.pseudoU.dist$rel_location, utr3.pseudoU.dist$rel_location)

metagene.cord <- c(m6a.metagene.coord, pseudoU.metagene.coord)
peakType <- c(rep("Nonsignificant peaks", length(m6a.metagene.coord)), 
         rep("Significant peaks", length(pseudoU.metagene.coord))) 
df <- data.frame(metagene.cord, peakType)

print("length of significant peaks")
print(length(pseudoU.metagene.coord))

print("length of nonsignificant peaks")
print(length(m6a.metagene.coord))

ggplot(df) + geom_density(aes(x = metagene.cord, colour = peakType), size = 1) + xlim(0, 3) + 
  theme_bw() + geom_vline(xintercept = 1:2, col = "grey") + ggtitle(fact) + theme(
    plot.title = element_text(color="red", size=24, face="bold.italic"))

ggsave(paste(fact, "metagene_plot.png", sep="_"), plot = last_plot(), device = "png", path = getwd(),
       scale = 1, width = 8, height = 8, units = "in",
       dpi = 300, limitsize = TRUE)

print(paste("finish metaPlotR for", fact))
