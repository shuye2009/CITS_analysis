#!/bin/bash
# This pipeline is used after CITS jobs have been finished 

scriptDir=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/scripts_program/CITS_CIMS

function combine(){
	
	local TARGET=$1 ## name of CLIP target or "Input" for Input 
	local OUTdir=$2 ## /home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/CLIP_cits  ## the location of the cits file
	local pv=$3 ## pvalue cutoff

	echo "Merge crosslink sites from replicates into one bed file"
	rm ${OUTdir}/*.*

	for citsd in $(ls -d CITS_${TARGET}*); do
                echo "processing ${citsd} for combined sites"
                factor=${citsd//CITS_/} #substiture CITS_ with an empty string
		${scriptDir}/CITS_analysis_thresholding.sh $factor ./$citsd $pv
	
		cat ./${citsd}/${factor}.tag.CITS.${pv}.bed >> ${OUTdir}/combined_CITS_${pv}_${TARGET}.bed
		cat ./${citsd}/${factor}.tag.CIMS.${pv}.bed >> ${OUTdir}/combined_CIMS_${pv}_${TARGET}.bed
                cat ./${citsd}/${factor}.crossLink_site.${pv}.bed >> ${OUTdir}/combined_crosslink_${pv}_${TARGET}.bed
        done


	echo "Filtering for unique crosslink sites using mergeBed"
	rn=$(echo $RANDOM | md5sum | head -c 20; echo)
	rd=/dev/shm/${TARGET}_$rn
	mkdir $rd

	cat ${OUTdir}/combined_CIMS_${pv}_${TARGET}.bed | awk -v OFS="\t" '{if($6!="+" && $6!="-")$6="+"; print $0}' | sort -T $rd -k1,1 -k2,2n | mergeBed -i stdin -s -d 0 -c 4,5,6 -o collapse,sum,distinct -delim "|" > ${OUTdir}/combined_CIMS_${TARGET}.uniq.bed
        cat ${OUTdir}/combined_CIMS_${TARGET}.uniq.bed | sort -k5,5nr > ${OUTdir}/combined_CIMS_${pv}_${TARGET}.uniq_sorted.bed
	
	cat ${OUTdir}/combined_CITS_${pv}_${TARGET}.bed | awk -v OFS="\t" '{if($6!="+" && $6!="-")$6="+"; print $0}' | sort -T $rd -k1,1 -k2,2n | mergeBed -i stdin -s -d 0 -c 4,5,6 -o collapse,sum,distinct -delim "|" > ${OUTdir}/combined_CITS_${TARGET}.uniq.bed
        cat ${OUTdir}/combined_CITS_${TARGET}.uniq.bed | sort -k5,5nr > ${OUTdir}/combined_CITS_${pv}_${TARGET}.uniq_sorted.bed

	cat ${OUTdir}/combined_crosslink_${pv}_${TARGET}.bed | awk -v OFS="\t" '{if($6!="+" && $6!="-")$6="+"; print $0}' | sort -T $rd -k1,1 -k2,2n | mergeBed -i stdin -s -d 0 -c 4,5,6 -o collapse,sum,distinct -delim "|" > ${OUTdir}/combined_crosslink_${TARGET}.uniq.bed
	cat ${OUTdir}/combined_crosslink_${TARGET}.uniq.bed | sort -k5,5nr > ${OUTdir}/combined_crosslink_${pv}_${TARGET}.uniq_sorted.bed

	rm -rf $rd
	rm ${OUTdir}/combined_CIMS_${TARGET}.uniq.bed ${OUTdir}/combined_CIMS_${pv}_${TARGET}.bed
	rm ${OUTdir}/combined_CITS_${TARGET}.uniq.bed ${OUTdir}/combined_CITS_${pv}_${TARGET}.bed
	rm ${OUTdir}/combined_crosslink_${TARGET}.uniq.bed ${OUTdir}/combined_crosslink_${pv}_${TARGET}.bed
}

function select_top(){	

	local TARGET=$1 ## name of CLIP target or "Input" for Input
        local OUTdir=$2 ## /home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/CLIP_cits  ## the location of the cits file
        local percent=$3 ## top percent to keep
	local pv=$4

	echo "Get top $percent percent of peaks"
	
	if (( "$percent" <= 100 )); then
		top=$(awk -v p=$percent 'END{printf("%.0f\n", FNR*p/100)}' ${OUTdir}/combined_CIMS_${pv}_${TARGET}.uniq_sorted.bed)
                head -n $top ${OUTdir}/combined_CIMS_${pv}_${TARGET}.uniq_sorted.bed > ${OUTdir}/combined_CIMS_${pv}_${TARGET}_top${percent}percent.bed

		top=$(awk -v p=$percent 'END{printf("%.0f\n", FNR*p/100)}' ${OUTdir}/combined_CITS_${pv}_${TARGET}.uniq_sorted.bed)
                head -n $top ${OUTdir}/combined_CITS_${pv}_${TARGET}.uniq_sorted.bed > ${OUTdir}/combined_CITS_${pv}_${TARGET}_top${percent}percent.bed

		top=$(awk -v p=$percent 'END{printf("%.0f\n", FNR*p/100)}' ${OUTdir}/combined_crosslink_${pv}_${TARGET}.uniq_sorted.bed)
		head -n $top ${OUTdir}/combined_crosslink_${pv}_${TARGET}.uniq_sorted.bed > ${OUTdir}/combined_crosslink_${pv}_${TARGET}_top${percent}percent.bed
	else
		top=$percent
		head -n $top ${OUTdir}/combined_CIMS_${pv}_${TARGET}.uniq_sorted.bed > ${OUTdir}/combined_CIMS_${pv}_${TARGET}_top${percent}.bed
		head -n $top ${OUTdir}/combined_CITS_${pv}_${TARGET}.uniq_sorted.bed > ${OUTdir}/combined_CITS_${pv}_${TARGET}_top${percent}.bed
		head -n $top ${OUTdir}/combined_crosslink_${pv}_${TARGET}.uniq_sorted.bed > ${OUTdir}/combined_crosslink_${pv}_${TARGET}_top${percent}.bed
	fi
}	

function main(){
        local TARGET=$1 ## name of CLIP target or "Input" for Input
        local TARGETdir=$2 ## /home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/CLIP_cits  ## the location of the cits file
        local pv=$3 ## pvalue cutoff
	local percent=$4

        if [ -z $pv ]; then pv=0.05; fi
	if [ -z $percent ]; then percent=10; fi

        cd ${TARGETdir}

        OUTdir=${TARGETdir}/${TARGET}_CITS_crosslinkSite
        mkdir -p $OUTdir

	if (( "$percent" <= 0 )) ; then
		## refresh combined_crosslink_${TARGET}.uniq_sorted.bed if it is present, otherwise make a new one
		combine "$TARGET" "$OUTdir" "$pv"
	elif [[ ! -s ${OUTdir}/combined_crosslink_${pv}_${TARGET}.uniq_sorted.bed ]]; then
		## test existing combined_crosslink_${TARGET}.uniq_sorted.bed, if it is empty, make a new one
		combine "$TARGET" "$OUTdir" "$pv"
		select_top "$TARGET" "$OUTdir" "$percent" "$pv"
	else
		select_top "$TARGET" "$OUTdir" "$percent" "$pv"
	fi
}

main "$@"
## END
