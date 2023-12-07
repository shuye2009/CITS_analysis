#!/bin/bash
# This pipeline is used after CITS jobs have been finished 
# Input should be processed first, the following procedures are performed:
# 1. combining Input crosslink sits into one bed file, 
# 2. filtering for unique crosslink sites, 
# 3. annotating crosslink sites relative to polyA motif, 
# 4. finding the genes where the crosslink sites are located, 
# 5. finding motif in all crossslink sites,
# 6. finding motif in 3UTR crosslink sites
#
# For a CLIP target, the following procedures are performed:
# 1. merge crosslink sites from replicates into one bed file, 
# 2. filtering for unique crosslink sites,
# 3. annotating crosslink sites relative to polyA motif,
# 4. finding the genes where the crosslink sites are located,
# 5. finding motif in all crossslink sites,
# 6. finding motif in 3UTR crosslink sites,
# 7. generating metagene plot
# 8. counting peaks in 3UTR

TARGET=$1 ## name of CLIP target or "Input" for Input 
TARGETdir=$2 ## /home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/CLIP_cits  ## the location of the cits file
TARGETbam=$3 ## /home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/CLIP_bam

INPUTdir=$4 ## /home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/Input_cits
INPUTbam=$5 ## /home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/Input_bam
scriptDir=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/scripts_program/CITS_CIMS
pv=1e-10
plotd=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/results_output/CITS_peak_in_gene_parts_${pv}
if test "$TARGET" = "Input"; then
	
	motifoutd=${INPUTdir}/homer_analysis
	mkdir -p $motifoutd

	if [ 1 -eq 1 ]; then
	### merge all input peaks into 1 bed file
	cd ${INPUTdir}
	echo "1. combining Input crosslink sits into one bed file"
	if [ -f combined_crosslink_Input.bed ]; then rm combined_crosslink_Input.bed; fi
	for inputd in $(ls -d CITS_Input*); do 
		echo "processing ${inputd}"
		factor=${inputd//CITS_/} #substiture CITS_ with an empty string
		${scriptDir}/CITS_analysis_thresholding.sh $factor ./$inputd $pv
		cat ./${inputd}/${factor}.crossLink_site.${pv}.bed >> combined_crosslink_Input.bed	
	done
	
	### collapse recurrent peaks into unique peaks
	echo "2. filtering for unique crosslink sites"
	cut -f1,2,3,6 combined_crosslink_Input.bed | sort | uniq -c | sort -k1,1nr > combined_crosslink_Input.bed.uniq.txt
	### output as a bed format file
	awk '{print $2"\t"$3"\t"$4"\tNA\t"$1"\t"$5}' combined_crosslink_Input.bed.uniq.txt > combined_crosslink_Input.uniq.bed
	fi

	if [ 1 -eq 1 ]; then
	echo "3. annotating crosslink sites relative to polyA motif"
	motifd=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/homer_analysis/motifs
        motif1=polyAsite.motif

	annotatePeaks.pl combined_crosslink_Input.uniq.bed hg19 -size 200 -norevopp -m ${motifd}/$motif1 -mdist -annStats Input_CITS_crosslinkSite_annotation_stats.txt > Input_CITS_crosslinkSite_${motif1}_annotation.txt

## homer annotatePeaks.pl only annotate genes of nearest TTS, not the genes where the peak is located. the perl program does this task. Also the crosslink sites located in 3UTRs are output in a separate bed file
        fi

	if [ 1 -eq 1 ]; then
	echo "4. finding the genes where the crosslink sites are located"
	perl_command=${scriptDir}/annotate_refseq_genes.pl
        perl ${perl_command} Input Input_CITS_crosslinkSite_${motif1}_annotation.txt ${INPUTdir}
	fi
	
	if [ 1 -eq 1 ]; then
## Find RNA motif around crosslink sites
	echo "5. finding motif in all crossslink sites"
        findmotif_command=${scriptDir}/homer_motif_from_CITS_sites_bc.sh
        ${findmotif_command} Input ${INPUTdir} ${motifoutd} ${INPUTdir} combined_crosslink_Input.uniq.bed

	echo "6. finding motif in 3UTR crosslink sites"
	findmotif_command=${scriptDir}/homer_motif_from_CITS_sites_in3UTR_bc.sh
        ${findmotif_command} Input ${INPUTdir} ${motifoutd} ${INPUTdir} ${TARGET}_CITS_crosslinkSite_in_3UTR.bed
	fi


else
	
	cd ${TARGETdir}

	OUTdir=${TARGETdir}/${TARGET}_CITS_crosslinkSite
	mkdir -p $OUTdir
	motifoutd=${OUTdir}/homer_analysis
        rm -r $motifoutd
if [ 1 -eq 1 ]; then
	echo "1. merge crosslink sites from replicates into one bed file"
	if [ -f ${OUTdir}/combined_crosslink_${TARGET}.bed ]; then rm ${OUTdir}/combined_crosslink_${TARGET}.bed; fi
        for citsd in $(ls -d CITS_${TARGET}*); do
                echo "processing ${sp1d} for combined sites"
                factor=${citsd//CITS_/} #substiture CITS_ with an empty string
		${scriptDir}/CITS_analysis_thresholding.sh $factor ./$citsd $pv
                cat ./${citsd}/${factor}.crossLink_site.${pv}.bed >> ${OUTdir}/combined_crosslink_${TARGET}.bed
        done

	echo "2. filtering for unique crosslink sites using mergeBed"
		cat ${OUTdir}/combined_crosslink_${TARGET}.bed | awk -v OFS="\t" '{if($6!="+" && $6!="-")$6="+"; print $0}' | sort -T ./ -k1,1 -k2,2n | mergeBed -i stdin -s -d 0 -c 4,5,6 -o collapse,sum,distinct -delim "|" > ${OUTdir}/combined_crosslink_${TARGET}.uniq.bed

		top10=$(awk 'END{printf("%.0f\n", FNR/10)}' ${OUTdir}/combined_crosslink_${TARGET}.uniq.bed)
		cat ${OUTdir}/combined_crosslink_${TARGET}.uniq.bed | sort -T ./ -k5,5nr | head -n $top10 > ${OUTdir}/combined_crosslink_${TARGET}.uniq_top10percent.bed

        #cut -f1,2,3,6 ${OUTdir}/combined_crosslink_${TARGET}.bed | sort | uniq -c | sort -k1,1nr > ${OUTdir}/combined_crosslink_${TARGET}.bed.uniq.txt

	#awk '{print $2"\t"$3"\t"$4"\tNA\t"$1"\t"$5}' ${OUTdir}/combined_crosslink_${TARGET}.bed.uniq.txt > ${OUTdir}/combined_crosslink_${TARGET}.uniq.bed
	#rm ${OUTdir}/combined_crosslink_${TARGET}.bed.uniq.txt
fi

if [ 1 -eq 0 ]; then
	echo "3. annotating peaks with respect to polyA motif using homer"
	motifd=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/homer_analysis/motifs
	motif1=polyAsite.motif
	annotatePeaks.pl ${OUTdir}/combined_crosslink_${TARGET}.uniq.bed hg19 -size 200 -norevopp -m ${motifd}/$motif1 -mdist -annStats ${OUTdir}/${TARGET}_CITS_crosslinkSite_annotation_stats.txt > ${OUTdir}/${TARGET}_CITS_crosslinkSite_${motif1}_annotation.txt
	
## homer annotatePeaks.pl only annotate genes of nearest TTS, not the genes where the peak is located. the perl program does this task.
	echo "4. annotating peaks with respect to genes"
	perl_command=${scriptDir}/annotate_refseq_genes.pl ## output ${factor}_CITS_crosslinkSite_in_3UTR.bed
	perl ${perl_command} ${TARGET} ${TARGET}_CITS_crosslinkSite_${motif1}_annotation.txt ${OUTdir}

	cut -f29,30 ${OUTdir}/${TARGET}_CITS_crosslinkSite_${motif1}_annotation_description.txt | grep "protein-coding" | cut -f1 | sort | uniq > ${OUTdir}/${TARGET}_CITS_crosslinkSite_${motif1}_annotation_description_topGO.txt
        cut -f29,30 ${OUTdir}/${TARGET}_CITS_crosslinkSite_${motif1}_annotation_description.txt | grep "ncRNA" | cut -f1 | sort | uniq > ${OUTdir}/${TARGET}_CITS_crosslinkSite_${motif1}_annotation_description_ncRNA.txt	
fi
if [ 1 -eq 0 ]; then
## Find RNA motif around crosslink sites
	echo "5. finding motif for crosslink sites"
	findmotif_command=${scriptDir}/homer_motif_from_CITS_sites_bc.sh
	## motifoutd=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/results_output/homer_analysis
	${findmotif_command} ${TARGET} ${OUTdir} ${motifoutd} ${INPUTdir} combined_crosslink_${TARGET}.uniq.bed

	echo "6. finding motif for crosslink sites in 3UTR"
        findmotif_command=${scriptDir}/homer_motif_from_CITS_sites_in3UTR_bc.sh
        ${findmotif_command} ${TARGET} ${OUTdir} ${motifoutd} ${INPUTdir} ${TARGET}_CITS_crosslinkSite_in_3UTR.bed
fi
if [ 1 -eq 0 ]; then  ## still applicable

	echo "7. producing metagene plot"
	plot_command=${scriptDir}/plot_CIMS_CITS_results_CLIP_generic.sh
	${plot_command} combined_crosslink_${TARGET}.uniq.bed ${OUTdir} combined_crosslink_Input.uniq.bed $INPUTdir 0,0,0,0,1,1 $plotd
fi
if [ 1 -eq 0 ]; then
	echo "8. count peaks in 3'UTR"
	bash ${scriptDir}/run_count_CITS_peaks_in_UTR3_PAU_siSP1.sh
fi
fi
echo "Finish all"

