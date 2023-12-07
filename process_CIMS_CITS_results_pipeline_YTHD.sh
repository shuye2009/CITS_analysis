#!/bin/bash
# This pipeline is used after CITS jobs have been finished 
# Input should be processed first, the following procedures are performed:
# 1. counting reads in Input bam files, 
# 2. combining Input crosslink sits into one bed file, 
# 3. filtering for unique crosslink sites, 
# 4. annotating crosslink sites relative to polyA motif, 
# 5. finding the genes where the crosslink sites are located, 
# 6. finding motif in all crossslink sites,
# 7. finding motif in 3UTR crosslink sites
#
# For a CLIP target, the following procedures are performed:
# 1. merge crosslink sites from replicates into one bed file, 
# 2. filtering for unique crosslink sites,
# 3. partition crosslink sites into two disjoint sets: one shared with Input, another one not shared with Input
# 4. perform DESeq2 analysis
# 5. annotating crosslink sites relative to polyA motif,
# 6. finding the genes where the crosslink sites are located,
# 7. finding motif in all crossslink sites,
# 8. finding motif in 3UTR crosslink sites,
# 9. computing distance matricies for metagene plot
# 10. generating metagene plot

TARGET=$1 ## name of CLIP target or "Input" for Input 
TARGETdir=$2 ## /home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/CLIP_cits  ## the location of the cits file
TARGETbam=$3 ## /home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/CLIP_bam

INPUTdir=$4 ## /home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/Input_cits
INPUTbam=$5 ## /home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/Input_bam
scriptDir=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/scripts_program/CITS_CIMS
pv=0.0001
plotd=/home/greenblattlab/shuyepu/Nabeel/YTHD/results_output/CITS_peak_in_gene_parts_${pv}
if test "$TARGET" = "Input"; then
	
	motifoutd=${INPUTdir}/homer_analysis
	mkdir -p $motifoutd

	cd ${INPUTbam}

	echo "1. counting reads in Input bam files"	
	### get input bam file size for use with DESeq2 analysis
	if [ -f Input_readCount_bam_files.tab ]; then rm Input_readCount_bam_files.tab; fi
        for f in Input*thUni.bam; do
                libsize=$(samtools view -c $f)
                echo -e "${f}\t${libsize}" >> Input_readCount_bam_files.tab
        done
	
	### merge all input peaks into 1 bed file
	cd ${INPUTdir}
	echo "2. combining Input crosslink sits into one bed file"
	if [ -f combined_crosslink_Input.bed ]; then rm combined_crosslink_Input.bed; fi
	for inputd in $(ls -d CITS_Input*); do 
		echo "processing ${inputd}"
		factor=${inputd//CITS_/} #substiture CITS_ with an empty string
		${scriptDir}/CITS_analysis_thresholding.sh $factor ./$inputd $pv
		cat ./${inputd}/${factor}.crossLink_site.${pv}.bed >> combined_crosslink_Input.bed	
	done

	### collapse recurrent peaks into unique peaks
	echo "3. filtering for unique crosslink sites"
	cut -f1,2,3,6 combined_crosslink_Input.bed | sort | uniq -c | sort -k1,1nr > combined_crosslink_Input.bed.uniq.txt
	### output as a bed format file
	awk '{print $2"\t"$3"\t"$4"\tNA\t"$1"\t"$5}' combined_crosslink_Input.bed.uniq.txt > combined_crosslink_Input.uniq.bed
	echo "4. annotating crosslink sites relative to polyA motif"
	motifd=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/homer_analysis/motifs
        motif1=polyAsite.motif

	annotatePeaks.pl combined_crosslink_Input.uniq.bed hg19 -size 200 -norevopp -m ${motifd}/$motif1 -mdist -annStats Input_CITS_crosslinkSite_annotation_stats.txt > Input_CITS_crosslinkSite_${motif1}_annotation.txt

## homer annotatePeaks.pl only annotate genes of nearest TTS, not the genes where the peak is located. the perl program does this task. Also the crosslink sites located in 3UTRs are output in a separate bed file
        
	echo "5. finding the genes where the crosslink sites are located"
	perl_command=${scriptDir}/annotate_refseq_genes.pl
        perl ${perl_command} Input Input_CITS_crosslinkSite_${motif1}_annotation.txt ${INPUTdir}
if [ 1 -eq 1 ]; then
## Find RNA motif around crosslink sites
	echo "6. finding motif in all crossslink sites"
        findmotif_command=${scriptDir}/homer_motif_from_CITS_sites_bc.sh
        ${findmotif_command} Input ${INPUTdir} ${motifoutd} ${INPUTdir} combined_crosslink_Input.uniq.bed

	echo "7. finding motif in 3UTR crosslink sites"
	findmotif_command=${scriptDir}/homer_motif_from_CITS_sites_in3UTR_bc.sh
        ${findmotif_command} Input ${INPUTdir} ${motifoutd} ${INPUTdir} ${TARGET}_CITS_crosslinkSite_in_3UTR.bed
fi


else
	
	cd ${TARGETdir}

	OUTdir=${TARGETdir}/${TARGET}_CITS_crosslinkSite
	mkdir -p $OUTdir
	motifoutd=${OUTdir}/homer_analysis
        mkdir -p $motifoutd
if [ 1 -eq 1 ]; then
	echo "1. merge crosslink sites from replicates into one bed file"
	if [ -f ${OUTdir}/combined_crosslink_${TARGET}.bed ]; then rm ${OUTdir}/combined_crosslink_${TARGET}.bed; fi
        for sp1d in $(ls -d CITS_${TARGET}*); do
                echo "processing ${sp1d} for combined sites"
                factor=${sp1d//CITS_/} #substiture CITS_ with an empty string
		#${scriptDir}/CITS_analysis_thresholding.sh $factor ./$sp1d $pv
                cat ./${sp1d}/${factor}_sorted.tag.CITS.${pv}.bed >> ${OUTdir}/tmp.bed
        done
	## remove chromosome "GL000*" and "MT" and add "chr" to chromosome column, output sorted bed based on tag count
	egrep -e "GL000|MT" -v ${OUTdir}/tmp.bed | awk '{print "chr"$0}' | sort -k5,5nr > ${OUTdir}/combined_crosslink_${TARGET}.bed
	rm -fr ${OUTdir}/tmp.bed
	## get top 25% peaks for motif discovery
	npeaks=$(wc -l ${OUTdir}/combined_crosslink_${TARGET}.bed | cut -d' ' -f1)
	fraction=$( echo $npeaks*0.25 | bc ) ## produce a float number
        n=$( echo ${fraction} | cut -d'.' -f1 )         ## convert the float number to integer
	head -n $n  ${OUTdir}/combined_crosslink_${TARGET}.bed | cut -f1,2,3,6 | sort | uniq -c | sort -k1,1nr > ${OUTdir}/combined_crosslink_${TARGET}.bed.top20pct.txt
	awk '{print $2"\t"$3"\t"$4"\tNA\t"$1"\t"$5}' ${OUTdir}/combined_crosslink_${TARGET}.bed.top20pct.txt > ${OUTdir}/combined_crosslink_${TARGET}.top20pct.uniq.bed
	rm ${OUTdir}/combined_crosslink_${TARGET}.bed.top20pct.txt

	echo "2. filtering for unique crosslink sites"
        cut -f1,2,3,6 ${OUTdir}/combined_crosslink_${TARGET}.bed | sort | uniq -c | sort -k1,1nr > ${OUTdir}/combined_crosslink_${TARGET}.bed.uniq.txt

	awk '{print $2"\t"$3"\t"$4"\tNA\t"$1"\t"$5}' ${OUTdir}/combined_crosslink_${TARGET}.bed.uniq.txt > ${OUTdir}/combined_crosslink_${TARGET}.uniq.bed
	rm ${OUTdir}/combined_crosslink_${TARGET}.bed.uniq.txt
fi

if [ 1 -eq 1 ]; then
	echo "5. annotating peaks with respect to polyA motif using homer"
	motifd=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/homer_analysis/motifs
	motif1=polyAsite.motif
	annotatePeaks.pl ${OUTdir}/combined_crosslink_${TARGET}.top20pct.uniq.bed hg19 -size 1 -strand + -norevopp -m ${motifd}/$motif1 -mdist -annStats ${OUTdir}/${TARGET}_CITS_crosslinkSite_annotation_stats.txt > ${OUTdir}/${TARGET}_CITS_crosslinkSite_${motif1}_annotation.txt
	
## homer annotatePeaks.pl only annotate genes of nearest TTS, not the genes where the peak is located. the perl program does this task.
	echo "6. annotating peaks with respect to genes"
	perl_command=${scriptDir}/annotate_refseq_genes.pl ## output ${factor}_CITS_crosslinkSite_in_3UTR.bed
	perl ${perl_command} ${TARGET} ${TARGET}_CITS_crosslinkSite_${motif1}_annotation.txt ${OUTdir}

	cut -f29,30 ${OUTdir}/${TARGET}_CITS_crosslinkSite_${motif1}_annotation_description.txt | grep "protein-coding" | cut -f1 | sort | uniq > ${OUTdir}/${TARGET}_CITS_crosslinkSite_${motif1}_annotation_description_topGO.txt
        cut -f29,30 ${OUTdir}/${TARGET}_CITS_crosslinkSite_${motif1}_annotation_description.txt | grep "ncRNA" | cut -f1 | sort | uniq > ${OUTdir}/${TARGET}_CITS_crosslinkSite_${motif1}_annotation_description_ncRNA.txt	
	## remove genes that are annotated as RNA genes, remove first line (header), output as bed format
	egrep -e "RNA" -v ${OUTdir}/${TARGET}_CITS_crosslinkSite_${motif1}_annotation.txt | tail -n +2 | awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$6"\t"$5}' > ${OUTdir}/combined_crosslink_${TARGET}.woRNA.uniq.bed

fi
if [ 1 -eq 1 ]; then
## Find RNA motif around crosslink sites
	echo "7. finding motif for crosslink sites"
	findmotif_command=${scriptDir}/homer_motif_from_CITS_sites_bc.sh
	## motifoutd=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/results_output/homer_analysis
	${findmotif_command} ${TARGET} ${OUTdir} ${motifoutd} ${INPUTdir} combined_crosslink_${TARGET}.woRNA.uniq.bed

	echo "8. finding motif for crosslink sites in 3UTR"
        findmotif_command=${scriptDir}/homer_motif_from_CITS_sites_in3UTR_bc.sh
        ${findmotif_command} ${TARGET} ${OUTdir} ${motifoutd} ${INPUTdir} ${TARGET}_CITS_crosslinkSite_in_3UTR.bed
fi
if [ 1 -eq 0 ]; then  ## still applicable

	echo "9. producing metagene plot"
	plot_command=${scriptDir}/plot_CIMS_CITS_results_CLIP_generic.sh
	${plot_command} combined_crosslink_${TARGET}.uniq.bed ${OUTdir} combined_crosslink_Input.uniq.bed $INPUTdir 1,1,0,0,0,0 $plotd
fi
if [ 1 -eq 0 ]; then
	echo "10. count peaks in 3'UTR"
	bash ${scriptDir}/run_count_CITS_peaks_in_UTR3_PAU_siSP1.sh
fi
fi
echo "Finish all"
