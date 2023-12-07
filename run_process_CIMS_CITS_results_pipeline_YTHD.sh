#! /bin/bash

cmd=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/scripts_program/CITS_CIMS/process_CIMS_CITS_results_pipeline_YTHD.sh
wd=/home/greenblattlab/shuyepu/Nabeel/YTHD/CLIP_cits
bamDir=/home/greenblattlab/shuyepu/Nabeel/YTHD/CLIP_bam
InputDir=/home/greenblattlab/shuyepu/Nabeel/YTHD/Input_cits
InputbamDir=/home/greenblattlab/shuyepu/Nabeel/YTHD/Input_bam
cd $wd

if [ 1 -eq 0 ]; then
	submitjob -m 40 $cmd Input $wd $bamDir $InputDir $InputbamDir
	sleep 5m
fi

if [ 1 -eq 1 ]; then
for f in $(ls -d CITS_YTHDF* | cut -d_ -f2 | cut -d. -f1 | sort | uniq | tr "\n" " "); do
	echo "submitting $f ..."	
	submitjob 100 -m 20 $cmd $f $wd $bamDir $InputDir $InputbamDir
	sleep 5s
	#$cmd $f $wd $bamDir $InputDir $InputbamDir
done 
fi
