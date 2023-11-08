#!/usr/bin/Rscript

## This script is used for filtering merged CITS peaks
## Only CITS with significantly more reads in CLIP files than in Input files are kept
## Author: Shuye Pu
## date created: June 17, 2022

print("Running $HOME/R_script/run_filter_peaks_by_Input.R ......")
source("/home/greenblattlab/shuyepu/Rscripts_on_diamond/Rscripts/GenomicPlot/R/GenomicPlot.R")

args <- commandArgs(trailingOnly = T)
queryfiles <- unlist(strsplit(args[1], split=",", fixed=T)) # CLIP files in bam format, separated by comma
inputfiles <- unlist(strsplit(args[2], split=",", fixed=T)) # Input files in bam format, separated by comma
peakfile <- args[3]   # a single peak file in bed format
inputdir <- args[4]
pcutoff <- args[5]
verbose <- args[6]

if(length(args) < 4) stop("query files, input files, peak file and input directory are mandatory!")
if(is.na(pcutoff)) pcutoff <- 0.5
if(is.na(verbose)) verbose <- FALSE

print(queryfiles)
print(inputfiles)
print(peakfile)
print(inputdir)
print(pcutoff)

results <- filter_peaks_by_differential(peakFile=peakfile, peakBam=queryfiles, refBam=inputfiles, refDir=inputdir, fix_width=21, genome="hg19", pcutoff=pcutoff) 

indices <- results$pos
peakbed <- read.delim(peakfile, header=F, sep="\t")
print(dim(peakbed))
peakbed <- peakbed[indices,]
print(dim(peakbed))
outfile <- gsub(".bed", "_filtered.bed", peakfile, fixed=T)
write.table(peakbed, outfile, sep="\t", col.names=F, row.names=F, quote=F)

if(verbose){
	da_table <- as.data.frame(results$res)
	outfile <- gsub(".bed", "_differential.tab", peakfile, fixed=T)
	write.table(da_table, outfile, sep="\t", col.names=F, row.names=F, quote=F)
}



## END
