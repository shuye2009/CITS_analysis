#!/bin/bash
## plot CITS peaks along with m6A peaks in common transcripts

scriptd=$HOME/Nabeel/clip_analysis/scripts_program/CITS_CIMS
m6as=( GSE79577_FTO_clusters.bed GSE102113_acRIP_peaks.bed M6Am_in_PCIF1_KD_lifted.bed m6A_CIMS_CITS.bed )
names=( "FTO" "AC4C" "m6Am" "m6A" )
m6ad=$HOME/Nabeel/ZBTB48_M6A/data
wd=$HOME/Nabeel/clip_analysis/data_input/CLIP_cits
cd $wd

for ((i=0; i<${#names[@]}; i++)); do
	m6a=${m6as[$i]}
	outd=$HOME/Nabeel/clip_analysis/results_output/CITS_peak_in_gene_parts_${names[$i]}
	mkdir -p $outd

	for d in *_CITS_crosslinkSite; do
		basef=$(basename $d _CITS_crosslinkSite)
		Rscript ${scriptd}/metaPlotR_CLIP_Input.R ${m6a},combined_crosslink_${basef}.uniq.bed ${m6ad},${wd}/${d} $outd TRUE
	done
done
