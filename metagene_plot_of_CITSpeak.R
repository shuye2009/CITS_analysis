## This script has to be called after "./run_deepTools_bamCompare_ratio_clip_pipeline.sh plot"
## has to be called after "deepTools_computeMatrix_metagene_parts.sh"
## For plotting single gene metagene profile, with or without a Input metagene profile

library(ComplexHeatmap)
library(circlize)

args = commandArgs(trailingOnly=TRUE)
targets <- args[1]
wds <- args[2] ## "/home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/SP1_analysis/mergedBW"
outd <- args[3]

targets <- unlist(strsplit(targets, ",", fixed=T))
wds <- unlist(strsplit(wds, ",", fixed=T))
regions <- c("prom", "utr5", "cds", "utr3", "tes")
region_len <- c(247, 46, 305, 155, 247) ## defined in "deepTools_computeMatrix_metagene_parts.sh"
region_target_mean_list <- list() ## a list of vectors

for(k in 1:length(wds)){
	target <- targets[k]
	wd <- wds[k]
	setwd(wd)
	region_target_mean <- NULL  ## a vector
	for(i in 1:length(regions)){
		region <- regions[i]
		len <- region_len[i]
		aFiles <- list.files(pattern=paste(target, "_over_", region,"_clip_matrix\\.gz$", sep=""))
		aFile <- aFiles[1]
		amatrix <- matrix(scan(aFile,what="character",skip=1),ncol=len+6, byrow=T)
		submatrix <- type.convert(amatrix[,7:ncol(amatrix)])
		submatrix[is.na(submatrix)] <- 0
		submatrix[is.nan(submatrix)] <- 0

		colm <- apply(submatrix, 2, mean) 
	
		region_target_mean <- c(region_target_mean, colm)
	}
	region_target_mean_list[[k]] <- region_target_mean

}

targetName <- paste(unlist(strsplit(targets[1], ".", fixed=T))[1:2], collapse="_")
if(grepl("combined", targets[1])){
	targetName <- unlist(strsplit(targets[1], "_", fixed=T))[3]
	targetName <- unlist(strsplit(targetName, ".", fixed=T))[1]
}

if(length(targets) == 2){
	targetName <- paste(targetName, "Input", sep="_")
}
xaxis <- seq(1, sum(region_len))
png(paste(outd, "/", targetName, "_peak_density_metagene.png", sep=""), width=8, height=8, unit="in", res=300)
yaxis <- scale(region_target_mean_list[[1]])
plot(smooth.spline(xaxis, yaxis, df=30), type="l", lwd=3, ylim=c(-2,2), col="cyan", xlab="Bins", ylab="Normalized CITS peak density", main=paste("CITS peaks of",targetName))

if(length(targets) == 2){
       	yaxis2 <- scale(region_target_mean_list[[2]])
	lines(smooth.spline(xaxis, yaxis2, df=30), lwd=3, col="orange")
	labels <- lapply(targets, function(x) unlist(strsplit(x, "_", fixed=T))[3])
	labels <- lapply(labels, function(x) unlist(strsplit(x, ".", fixed=T))[1])
	legend("topright", legend=labels, col=c("cyan", "orange"), lwd=3)
}

mtext(text=c("TSS-1kb", "5'UTR", "CDS", "3'UTR", "TES+2kb"), side=3, adj=c(0.1, 0.27, 0.45, 0.67, 0.9))
	#mtext(text=c("5'UTR", "CDS", "3'UTR"), side=3, adj=c(0.02, 0.34, 0.8))
abline(v=c(247, 293, 598, 753), lty=2)
dev.off()
  

## plot intron separately, currently not used, need update if put to use
if(0){
regions <- c("intron")
region_len <- c(1000)
region_target_mean <- NULL ## a dataframe
for(i in 1:length(regions)){
        region <- regions[i]
        len <- region_len[i]
        inputFiles <- list.files(pattern=paste(".*", norm, "_over_", region,"_clip_matrix\\.gz$", sep=""))

        target_mean <- list()
        for(aFile in inputFiles){
                #aFile <- inputFiles[1]
                target <- unlist(strsplit(aFile, ".", fixed=T))[1]
                amatrix <- matrix(scan(aFile,what="character",skip=1),ncol=len+6, byrow=T)
                submatrix <- type.convert(amatrix[,7:ncol(amatrix)])
                submatrix[is.na(submatrix)] <- 0
                submatrix[is.nan(submatrix)] <- 0

                colm <- apply(submatrix, 2, mean)
                target_mean[[target]] <- colm
        }
        target_mean_df <- as.data.frame(target_mean)  ## target as col names
        region_target_mean <- rbind(region_target_mean, target_mean_df)

}

region_target_matrix <- t(scale(as.matrix(region_target_mean)))
region_target_matrix[is.na(region_target_matrix)] <- 0
region_target_matrix[is.nan(region_target_matrix)] <- 0
lowb <- min(region_target_matrix)
highb <- max(region_target_matrix)
midb <- (lowb+highb)/2

type <- rep(regions, region_len)
ha = HeatmapAnnotation(df = data.frame(feature = type))

png(paste("CLIP_signal_intron", norm, "heatmap.png", sep="_"))
print(Heatmap(region_target_matrix, name="NormalizedLog2Rato", col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), top_annotation = ha, show_row_names = TRUE, show_column_names = FALSE, cluster_columns=FALSE, row_names_side="left"))
dev.off()
}



