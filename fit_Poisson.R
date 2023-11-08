#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = T)

countFile <- args[1]
pcutoff <- as.numeric(args[2])
prefix <- args[3]

print(paste("countFile", countFile))
print(paste("pcutoff", pcutoff))
print(paste("prefix", prefix))

count_table <- as.data.frame(data.table::fread(countFile))
print(paste(c("dim count_table", dim(count_table))))

print(summary(count_table))

if(pcutoff < 1){
	counts <- as.numeric(count_table[,4])
	counts <- counts[!is.na(counts)]
	uniq_counts <- sort(unique(counts))
	print(paste0("fitting Poisson model for ", prefix))
	para_Poisson <- MASS::fitdistr(counts, "Poisson")
	#para_NB <- MASS::fitdistr(counts, "negative binomial")

	dPoisson <- dpois(x=uniq_counts, lambda=para_Poisson$estimate, log=F)
	dt <- data.frame(uniq_counts, dPoisson)

	print(paste0("finding cutoff for probability ", pcutoff))
	cv <- qpois(pcutoff, lambda=para_Poisson$estimate, lower.tail=F, log.p=F)
}else{
	cv <- pcutoff
}

hotspot <- count_table[as.numeric(count_table[,4]) > cv,]

id_col <- paste0("hotspot", seq.int(nrow(hotspot)))
hotspot <- cbind(hotspot[,1:3], id_col, hotspot[,4:5])

print(paste(c("dim hotspot", dim(hotspot))))
dir <- dirname(countFile)

## output crosslink hotspots such that p(count > cv) < pcutoff
write.table(hotspot, file.path(dir,paste0(prefix,"_crosslink_hotspot.bed")), col.names=F, row.names=F, sep="\t", quote=F)

## output probability density of observed counts based on estimated lambda
if(pcutoff < 1) write.table(dt,  file.path(dir, paste0(prefix,"_Poisson_fitted_pvalue.tab")), col.names=T, row.names=F, sep="\t", quote=F)

print(paste("The cutoff for crosslink hotspot is:", cv)) 
