#!/bin/bash

factor=$1
peakbed=$2
peakdir=$3
outd=$4
pau=$5 ## 0_100 or 25_75
scriptd=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/scripts_program/CITS_CIMS
utrd=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/results_output/GO_results_${pau}

mkdir -p $outd
#peakbed=combined_crosslink_${factor}.uniq.bed
changes=(Up NoChange Down)
if [ -f $outd/${factor}_peaks_count_table_${pau}.tab ]; then rm $outd/${factor}_peaks_count_table_${pau}.tab; fi

for utr in putr dutr; do
	utrs=(PAU_analysis_siSP1_up_list_${utr}.bed PAU_analysis_siSP1_noChange_list_${utr}.bed PAU_analysis_siSP1_down_list_${utr}.bed)

	for((i=0; i<${#changes[*]}; i++)); do

		change=${changes[i]}
		utrbed=${utrs[i]}
		bedtools intersect -b $peakdir/$peakbed -a $utrd/$utrbed -s -c | awk -v autr="$utr" -v achange="$change" 'BEGIN{OFS="\t"}{print $0,achange,autr}' >> $outd/${factor}_peaks_count_table_${pau}.tab
	
	# shrink utr to the last 1000 nt of the 3', not working because the peak density is too low
	#	awk 'BEGIN{OFS="\t"}{if($6=="+"){print $1,$3-1000,$3,$4,$5,$6,$7,$8}else{print $1,$2,$2+1000,$4,$5,$6,$7,$8}}' $utrd/$utrbed | bedtools intersect -b $peakdir/$peakbed -a stdin -s -c | awk -v autr="$utr" -v achange="$change" 'BEGIN{OFS="\t"}{print $0,achange,autr}' > $outd/${factor}_peaks_count_table1000_${pau}.tab
	done
done

## END
