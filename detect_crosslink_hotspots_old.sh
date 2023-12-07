#!/bin/bash

function collect_CITS(){
	echo "merging files"
	local typeCL=$1
	local i=0
	for d in $(ls -d CITS_*); do 
		cd $d 
		fa=${d//CITS_/} 

		#rm ${fa}.tag.bed ## these files are an intermediates for CIMS and CITS peak calling, once done, they are never used afterwards
		#rm ${fa}.mutation.txt
		#rm ${fa}.tag.mutation.txt
		if [[ -e ${fa}.tag.del.bed ]]; then rm ${fa}.tag.del.bed ${fa}.tag.sub.bed ${fa}.tag.ins.bed; fi
		if [[ ! $fa =~ .*m6A.* ]]; then
			if [[ ! -s "${fa}.tag.CIMS.bed" ]]; then # -s means file not empty
				cat  ${fa}.tag.sub.CIMS.bed ${fa}.tag.del.CIMS.bed > ${fa}.tag.CIMS.bed
                        fi

			if [[ -s "${fa}.tag.${typeCL}.bed" ]]; then 
				echo "$d $i"
				cat ${fa}.tag.${typeCL}.bed >> ../Hotspot/all_factor.${typeCL}.bed 
				i=$((i+1))
			fi
		fi
		cd .. 
	done
	all=$i
}

function count_peak_frequency(){
	echo "counting peak frequency"
	local typeCL=$1
	cat ./Hotspot/all_factor.${typeCL}.bed | sort -T /dev/shm/ -k1,1 -k2,2n | mergeBed -i stdin -s -d -1 -c 4,5,6 -o collapse,count,distinct -delim "|" | cut -f1-3,5-7 > ./Hotspot/all_factor.${typeCL}_count.bed	
	rm ./Hotspot/all_factor.${typeCL}.bed
}

function main(){
	Rcmd=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/scripts_program/CITS_CIMS/fit_Poisson.R
	wd=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/CLIP_cits
	cd $wd
	mkdir -p $wd/Hotspot
	typeCL=$1 ## CITS or CIMS
	all=0

	if [[ -e $wd/Hotspot/all_factor.${typeCL}.bed ]]; then rm $wd/Hotspot/all_factor.${typeCL}.bed; fi

	collect_CITS "$typeCL"
	half=$(($all/2 + 1))
	echo $all $half

	count_peak_frequency "$typeCL"

	Rscript $Rcmd $wd/Hotspot/all_factor.${typeCL}_count.bed $half $typeCL
}

main "$@"
