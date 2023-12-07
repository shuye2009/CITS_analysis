#!/bin/bash

cmd=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/scripts_program/CITS_CIMS/process_CIMS_CITS_results_pipeline_generic.sh

function main(){
	local wd=$1 ##/home/greenblattlab/shuyepu/Nabeel/CLIP2021/CLIP_cits
	local pv=$2
	local factor=$3

	if [[ $factor == "all" ]]; then $factor=""; fi

	if [[ ! $# -eq 3 ]]; then 
		echo "Incorrect number of arguments: <working_directory> <pvalue> <factor>"
		exit 1
	fi

	cd $wd

	for f in $(ls -d CITS_${factor}* | cut -d_ -f2 | cut -d. -f1 | cut -d_ -f1 | sort | uniq | tr "\n" " "); do
		echo "processing $f ..."
		for top in 0 10 20; do	
			$cmd $f $wd $pv $top
		done
	done 
}

wd=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/CLIP_cits
pv=0.001
factor=$1
if [[ ! $# -eq 1 ]]; then
                echo "Incorrect number of arguments: <factor>"
                exit 1
        fi

main "$wd" "$pv" "$factor"


