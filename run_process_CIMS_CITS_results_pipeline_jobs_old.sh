#! /bin/bash

cmd=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/scripts_program/CITS_CIMS/process_CIMS_CITS_results_pipeline_generic.sh
wd=/home/greenblattlab/shuyepu/Nabeel/CLIP2021/CLIP_cits
bamDir=/home/greenblattlab/shuyepu/Nabeel/CLIP2021/CLIP_bam
InputDir=/home/greenblattlab/shuyepu/Nabeel/CLIP2021/Input_cits
InputbamDir=/home/greenblattlab/shuyepu/Nabeel/CLIP2021/Input_bam
cd $wd

if [ 1 -eq 0 ]; then
	submitjob -m 20 $cmd Input $wd $bamDir $InputDir $InputbamDir
	sleep 5m
fi

if [ 1 -eq 1 ]; then
for f in $(ls -d CITS_*m6A* | cut -d_ -f2 | cut -d. -f1 | cut -d_ -f1 | sort | uniq | tr "\n" " "); do
	echo "submitting $f ..."	
	#submitjob 10 -m 20 $cmd $f $wd $bamDir $InputDir $InputbamDir
	#sleep 5s
	$cmd $f $wd $bamDir $InputDir $InputbamDir
done 
fi
