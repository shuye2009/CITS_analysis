#! /bin/bash

cmd=$HOME/Nabeel/clip_analysis/scripts_program/CITS_CIMS/count_CITS_peaks_in_UTR3_PAU_siSP1.sh

if [ 1 -eq 0 ]; then
dir=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/SP1_analysis_AB/CLIP_cits
outd=$HOME/Nabeel/clip_analysis/results_output/SP1_analysis_AB/PAU_siSP1_peak_count_in_3UTR_PAU
cd $dir

for d in $(ls -d *_crosslinkSite); do
	f=$(echo "$d" | cut -d_ -f1)
	echo "$d $f"
	$cmd $f combined_crosslink_${f}.uniq.bed $dir/$d $outd
done 
fi

if [ 1 -eq 0 ]; then
dir=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/SP1_analysis_AB/Input_cits
outd=$HOME/Nabeel/clip_analysis/results_output/SP1_analysis_AB/PAU_siSP1_peak_count_in_3UTR_PAU
cd $dir
	$cmd Input combined_crosslink_Input.uniq.bed $dir $outd
fi
## END

if [ 1 -eq 1 ]; then
Rcmd1=$HOME/Nabeel/clip_analysis/scripts_program/CITS_CIMS/peak_targeted_UTR.R
Rcmd2=$HOME/Nabeel/clip_analysis/scripts_program/CITS_CIMS/plot_peak_count_in_UTR.R
Rcmd3=$HOME/Nabeel/clip_analysis/scripts_program/CITS_CIMS/plot_motif_frequency.R
dreme=$HOME/Nabeel/clip_analysis/data_input/SP1_2021/pureclip_peaks/dreme_SP1_S2S3Ecombined_crosslinkregions_wInputwCL/dreme.txt
dir=$HOME/Nabeel/clip_analysis/data_input/SP1_2021/pureclip_peaks
wd=$HOME/Nabeel/clip_analysis/data_input/SP1_2021/motif_frequency
plotd=$HOME/Nabeel/clip_analysis/data_input/SP1_2021/plots/motif_frequency
target=SP1_S2S3E
motif=DGGDRG
cd $wd
	for range in 0_100 25_75; do
		$cmd ${target} SP1_S2S3E_crosslinkregions_wInputwCL.bed $dir $wd $range
		Rscript $Rcmd1 $range SP1_S2S3E_peaks_count_table_${range}.tab $wd ${target}
		Rscript $Rcmd2 SP1_S2S3E_peaks_count_table_${range}.tab $wd ${target} $range $plotd

		cp SP1_S2S3E_peaks_count_table_${range}.tab SP1_S2S3E_peaks_count_table_${range}.bed
		count_motif_frequency_for_bed_bc.sh SP1_S2S3E_peaks_count_table_${range}.bed $motif $dreme false
	done

	for utr in dutr putr; do
		for pau in Down NoChange Up; do
			for range in 0_100 25_75; do
				for targetType in targeted nonTargeted; do
	
					count_motif_frequency_for_bed_bc.sh siSP1_PAU_${pau}_${utr}_${range}_${targetType}_by_${target}.bed $motif $dreme false
				done
			done
		done
	done

	Rscript $Rcmd3 $wd $motif $target $plotd
fi

## END
