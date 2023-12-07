#!/bin/bash

cmd=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/scripts_program/CITS_CIMS/process_CIMS_CITS_results_pipeline_generic.sh

function main(){
	local wd=$1 ##/home/greenblattlab/shuyepu/Nabeel/CLIP2021/CLIP_cits
	local pv=$2
	local factor=$3
	local fun=$4
	local replace=$5
	
	if [[ -z $replace ]]; then replace=false; fi
	if [[ $factor == "all" ]]; then factor="*"; fi

	if [[ $# -lt 4 ]]; then 
		echo "Incorrect number of arguments: <working_directory> <pvalue> <factor> <subfunction>"
		exit 1
	fi

	cd $wd

	shopt -s extglob ## to allow using of regex in ls
	for f in $(ls -d CITS_${factor}+(_|.)* | cut -d_ -f2 | cut -d. -f1 | cut -d_ -f1 | sort | uniq | tr "\n" " "); do
		if [[ $factor == "*" ]] && [[ $f != "Input" ]]; then
		 	if [[ ! -d ${f}_CITS_crosslinkSite ]] || [[ $replace = true ]]; then
				echo "processing $f ..."
				submitjob -w 48 -m 40 $cmd $f $wd $pv $fun $replace
			else
				echo "${f}_CITS_crosslinkSite exists already, to modify files in the directory, set 'replace'=true"
			fi
		elif [[ $factor != "*" ]]; then 
			echo "processing $f ..."
                                submitjob -w 100 -m 100 $cmd $f $wd $pv $fun $replace
		fi
	done 
}

wd=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/CLIP_cits

factor=$1
fun=$2
replace=$3 ## replace existing files
pv=$4 ## something like 10e-10 or 0.01

#if [[ $factor == "Input" ]]; then pv=10e-10; fi
if [[ $# -lt 2 ]]; then
        echo "Enter: <factor> for the first argument, use any factor but not 'all' for function crosslink_hotspot!"
		echo "Enter: combine, crosslink_hotspot, filter_peaks or discover_motif for the second argument"
		echo "Enter: true or false for the third argument, to indicate whether exiting files should be replaced"
		echo "Enter: pvalue to peak thresholding, like 0.01 or 10e-10"
        exit 1
fi

main "$wd" "$pv" "$factor" $fun $replace


