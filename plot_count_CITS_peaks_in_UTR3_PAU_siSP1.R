# normalize the count by length and plot a boxplot

args = commandArgs(trailingOnly=TRUE)

inputFiles <- args[1]
outputFile <- args[2]
wd <- args[3] ## "/home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/SP1_analysis/mergedBW"
labeles <- args[4]
target <- args[5]
#setwd("c:/RSYNC/Nabeel/ZBTB48_M6A/results/metaPlotR_plot")
#normalizePath(wd, winslash = '\\')
setwd(wd)

infiles <- unlist(strsplit(inputFiles, ",", fixed=T))
lab <- unlist(strsplit(labeles, ",", fixed=T))

count_list <- list()

for(i in 1:length(infiles)){

        my.dist.table <- type.convert(read.delim(infiles[i], header = F, stringsAsFactors = F))
        norm_count <- 1000*my.dist.table[,9]/(my.dist.table[,3] - my.dist.table[,2]) ## counts per kb
        count_list[[lab[i]]] <- norm_count
}

if(grepl("pdf", outputFile)){
	pdf(paste(wd,outputFile,sep="/"), width=8, height=6)
}else{
	png(paste(wd,outputFile,sep="/"), width=8, height=6, unit="in", res=300)
}

boxplot(count_list, ylab="peaks per kilobase", outline=F, main=paste(target, ": t-test p-values", sep=""))
ts <- t.test(count_list[["UP"]], count_list[["DOWN"]])
pv <- sprintf("%.3f", ts$p.value)
mtext(paste("UP vs. DOWN ", pv),  adj = 0, col = "red" )
ts <- t.test(count_list[["NOCHANGE"]], count_list[["DOWN"]])
pv <- sprintf("%.3f", ts$p.value)
mtext(paste("NOCHANGE vs. DOWN ", pv),  adj = 0.5, col = "red" )
ts <- t.test(count_list[["UP"]], count_list[["NOCHANGE"]])
pv <- sprintf("%.3f", ts$p.value)
mtext(paste("UP vs. NOCHANGE ", pv),  adj = 1.0, col = "red" )

dev.off()

if(0){
density_file <- gsub(".pdf", "_density.pdf", outputFile, fixed=T)
pdf(paste(wd,density_file,sep="/"), width=8, height=6)

plot(density(count_list[["NOCHANGE"]]), log="x", lwd=3, main=target, col="orange")
lines(density(count_list[["DOWN"]]), log="x", lwd=3, col="blue")
lines(density(count_list[["UP"]]), log="x", lwd=3, col="cyan")
legend("topright", legend=c("noChange", "Down", "Up"), lwd=3, col=c("orange", "blue", "cyan"))
dev.off()


if(grepl("pdf", outputFile)){
	quantile_file <- gsub(".pdf", "_quantile.pdf", outputFile, fixed=T)
	pdf(paste(wd,quantile_file,sep="/"), width=8, height=6)
}else{
	quantile_file <- gsub(".png", "_quantile.png", outputFile, fixed=T)
        png(paste(wd,quantile_file,sep="/"), width=8, height=6, unit="in", res=300)
}

probs <- seq(0, 1, 0.01)
plot(quantile(count_list[["NOCHANGE"]], probs=probs), probs, type="l", xlab="peaks per kilobase", ylab="quantile", lwd=3, main=target, col="orange")
lines(quantile(count_list[["DOWN"]], probs=probs), probs, lwd=3, col="blue")
lines(quantile(count_list[["UP"]], probs=probs), probs, lwd=3, col="cyan")
legend("bottomright", legend=c("noChange", "Down", "Up"), lwd=3, col=c("orange", "blue", "cyan"))


ks <- ks.test(quantile(count_list[["UP"]], probs=probs), quantile(count_list[["DOWN"]], probs=probs))
pv <- sprintf("%.3f", ks$p.value)
mtext(paste("UP vs. DOWN p-value = ",pv, " (ks-test)", sep=""),  adj = 0.5, col = "red" )

dev.off()
}
