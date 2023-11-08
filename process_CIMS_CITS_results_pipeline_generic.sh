#!/bin/bash
# This pipeline is used after CITS jobs have been finished 
# combine function collect CITS from individual samples into combined file, then merge overlapping sites into premerged.bed
# crosslink_hotspot function identify hotspots by examine all premerged bed files, CITS that appeary in more than 50% TARGETs are defined as hotspots
# filter_peaks function first remove peaks that are in hotspots, producing merged.bed, then remove peaks that have pvalue > 0.5 in differential analysis against Input at read level, producing merger_filtered.bed
# discover_motif function identify hexamer motifs using DREME

scriptDir=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/scripts_program/CITS_CIMS
bamDir=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/CLIP_bam
InputDir=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/Input_bam
Inputbam=$(ls $InputDir/*.bam | tr "\n" ",")

function combine(){
	
	local TARGET=$1 ## name of CLIP target or "Input" for Input 
	local TARGETdir=$2 ## /home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/CLIP_cits  ## the location of the cits file
	local pv=$3 ## pvalue cutoff
	local replace=$4

	cd $TARGETdir

	OUTdir=${TARGETdir}/${TARGET}_CITS_crosslinkSite
	mkdir -p $OUTdir

	echo "Merge crosslink sites from replicates into one bed file"
	if [[ $replace ]]; then rm -rf ${OUTdir}/*.*; fi

	rn=$(echo $RANDOM | md5sum | head -c 20; echo)
        rd=/dev/shm/${TARGET}_$rn
        mkdir $rd

	shopt -s extglob  # to allow pattern match in ls command
	for citsd in $(ls -d CITS_${TARGET}+(_|.)*); do
                echo "processing ${citsd} for combined sites"
                factor=${citsd//CITS_/} #substiture CITS_ with an empty string
		d=${TARGETdir}/$citsd
		if [[ -e ${d}/${factor}.tag.del.bed ]]; then rm ${d}/${factor}.tag.del.bed ${d}/${factor}.tag.sub.bed ${d}/${factor}.tag.ins.bed; fi

		if [[ ! -s "${d}/${factor}.tag.CIMS.bed" ]]; then # -s means file not empty
                	cat  ${d}/${factor}.tag.sub.CIMS.bed ${d}/${factor}.tag.del.CIMS.bed > ${d}/${factor}.tag.CIMS.bed
                fi

		${scriptDir}/CITS_analysis_thresholding.sh $factor $d $pv

		cat ${d}/${factor}.tag.CITS.${pv}.bed | sort -T $rd -k1,1 -k2,2n | mergeBed -i stdin -s -d 1 -c 4,5,6 -o collapse,count,distinct -delim "|" | cut -f1-6 >> ${OUTdir}/combined_CITS_${pv}_${TARGET}.bed
             	cat ${d}/${factor}.tag.CIMS.${pv}.bed | sort -T $rd -k1,1 -k2,2n | mergeBed -i stdin -s -d 1 -c 4,5,6 -o collapse,count,distinct -delim "|" | cut -f1-6 >> ${OUTdir}/combined_CIMS_${pv}_${TARGET}.bed

	done


	echo "Merging crosslink sites using mergeBed"

	cat ${OUTdir}/combined_CIMS_${pv}_${TARGET}.bed | sort -k1,1 -k2,2n | mergeBed -i stdin -s -d 1 -c 4,5,6 -o collapse,count,distinct -delim "&" | cut -f1-6 | sort -k5,5nr > ${OUTdir}/combined_CIMS_${pv}_${TARGET}.premerged.bed

	cat ${OUTdir}/combined_CITS_${pv}_${TARGET}.bed | sort -k1,1 -k2,2n | mergeBed -i stdin -s -d 1 -c 4,5,6 -o collapse,count,distinct -delim "&" | cut -f1-6 | sort -k5,5nr > ${OUTdir}/combined_CITS_${pv}_${TARGET}.premerged.bed

	cat ${OUTdir}/combined_CIMS_${pv}_${TARGET}.bed ${OUTdir}/combined_CITS_${pv}_${TARGET}.bed | sort -k1,1 -k2,2n | mergeBed -i stdin -s -d 1 -c 4,5,6 -o collapse,count,distinct -delim "&" | cut -f1-6 | sort -k5,5nr > ${OUTdir}/combined_crosslink_${pv}_${TARGET}.premerged.bed

	rm -rf $rd
}


function crosslink_hotspot(){

        Rcmd=$scriptDir/fit_Poisson.R
        local target=$1 ## can be anything, will be replaced by *
        local wd=$2 ## /home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/CLIP_cits  ## the location of the cits file
        local pv=$3 ## pvalue cutoff
	local replace=$4

        cd $wd
        mkdir -p $wd/Hotspot
        shopt -s extglob  # to allow pattern match in ls command

        for typeCL in CITS CIMS crosslink; do
                echo "merging files $typeCL"

                local i=0

                for TARGET in $(ls -d CITS_*+(_|.)* | cut -d_ -f2 | cut -d. -f1 | cut -d_ -f1 | sort | uniq | tr "\n" " "); do
                        if [[ ! $TARGET =~ .*m6A.* ]]; then

                                d=${TARGET}_CITS_crosslinkSite

                                if [[ -s "${d}/combined_${typeCL}_${pv}_${TARGET}.premerged.bed" ]]; then
                                        echo "$d $i"
                                        cat ${d}/combined_${typeCL}_${pv}_${TARGET}.premerged.bed >> $wd/Hotspot/all_factor.${typeCL}.bed
                                        i=$((i+1))
                                fi
                        fi
                done

                all=$i

                echo "counting peak frequency"
		cat $wd/Hotspot/all_factor.${typeCL}.bed | sort -T /dev/shm/ -k1,1 -k2,2n | mergeBed -i stdin -s -d -1 -c 4,5,6 -o collapse,count,distinct -delim "|" | cut -f1-3,5-7 > $wd/Hotspot/all_factor.${typeCL}_count.bed
                rm $wd/Hotspot/all_factor.${typeCL}.bed

                echo "defining hotspots"

                half=$(($all/2 + 1))
                echo $all $half

                Rscript $Rcmd $wd/Hotspot/all_factor.${typeCL}_count.bed $half $typeCL
        done
}



function filter_peaks(){

        local TARGET=$1 ## name of CLIP target or "Input" for Input
        local TARGETdir=$2 ## /home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/CLIP_cits  ## the location of the cits file
        local pv=$3 ## pvalue cutoff
	local replace=$4

        OUTdir=${TARGETdir}/${TARGET}_CITS_crosslinkSite
        mkdir -p $OUTdir

	echo "removing hotspots"

	rn=$(echo $RANDOM | md5sum | head -c 20; echo)
        rd=/dev/shm/${TARGET}_$rn
        mkdir $rd

	for typeCL in CITS CIMS crosslink; do
		hotspot=${TARGETdir}/Hotspot/${typeCL}_crosslink_hotspot.bed

		cat ${OUTdir}/combined_${typeCL}_${pv}_${TARGET}.premerged.bed | subtractBed -A -s -a stdin -b ${hotspot} | sort -T $rd -k1,1 -k2,2n | cut -f1-6 > ${OUTdir}/combined_${typeCL}_${pv}_${TARGET}.merged.bed
		cat ${OUTdir}/combined_${typeCL}_${pv}_${TARGET}.merged.bed | awk -v OFS="\t" '{if($5 > 1) print $0}' > ${OUTdir}/combined_${typeCL}_${pv}_${TARGET}.recurring.bed
	done

	rm -rf $rd

	echo "Collecting CLIP bam files"

        bamfiles=""
        shopt -s extglob  # to allow pattern match in ls command
        for citsd in $(ls -d CITS_${TARGET}+(_|.)*); do
                echo "processing ${citsd} for bam file name"
                factor=${citsd//CITS_/} #substiture CITS_ with an empty string

                if [[ -z $bamfiles ]]; then
                        bamfiles="$bamDir/${factor}.bam"
                else
                        bamfiles="$bamfiles,$bamDir/${factor}.bam"
                fi

	done
 	if [[ $TARGET == "Input" ]]; then bamfiles=$Inputbam; fi
	echo "Filtering merged peaks using input read counts"
        
	Rscript $HOME/R_script/run_filter_peaks_by_Input.R $bamfiles $Inputbam ${OUTdir}/combined_CITS_${pv}_${TARGET}.merged.bed $InputDir 0.5 FALSE

	rm ${OUTdir}/combined_CITS_${pv}_${TARGET}.bed ${OUTdir}/combined_CIMS_${pv}_${TARGET}.bed
}

function discover_motif(){

        local TARGET=$1 ## name of CLIP target or "Input" for Input
        local TARGETdir=$2 ## /home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/CLIP_cits  ## the location of the cits file
        local pv=$3 ## pvalue cutoff
	local replace=$4

        OUTdir=${TARGETdir}/${TARGET}_CITS_crosslinkSite
	cd $OUTdir
        echo "Discovering motif with DREME"
     
	## bed intervals will be extended 10 bp on each side
	#crosslinksites_to_dreme_motif.sh combined_CIMS_${pv}_${TARGET}.recurring.bed 10
	crosslinksites_to_dreme_motif.sh combined_CITS_${pv}_${TARGET}.merged_filtered.bed 10
	#crosslinksites_to_dreme_motif.sh combined_crosslink_${pv}_${TARGET}.recurring.bed 10
	
}


function main(){
        local TARGET=$1 ## name of CLIP target or "Input" for Input
        local TARGETdir=$2 ## /home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/CLIP_cits  ## the location of the cits file
        local pv=$3 ## pvalue cutoff
	local fun=$4 ## name of function to call
	local replace=$5
        if [ -z $pv ]; then pv=0.05; fi

        cd ${TARGETdir}

	$fun "$TARGET" "$TARGETdir" "$pv" $replace
}

main "$@"
## END
