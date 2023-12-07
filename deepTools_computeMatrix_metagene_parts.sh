#!/bin/bash
# compute metagene matrix for UTR, CDS, TES and introns. 
# plotting will be done in heatmap_of_clip_signal_in_parts.R

f=$1 # .bw file
wd=$2 # directory of .bw file
outd=$3 ## /home/greenblattlab/shuyepu/Nabeel/clip_analysis/results_output/clip_signal_in_gene_parts/
#mkdir -p ${outd}

DIR=$HOME/genomic_feature
prom=${DIR}/gencode.v19.protein_coding_transcript_promoter.bed
utr3=${DIR}/gencode.v19.protein_coding_transcript_UTR3.bed12
utr5=${DIR}/gencode.v19.protein_coding_transcript_UTR5.bed12
cds=${DIR}/gencode.v19.protein_coding_transcript_CDS.bed12
tes=${DIR}/non_overlapping_protein_coding_genes_eag_500-500_500_0-2000_hg19.bed
intron=${DIR}/gencode.v19.protein_coding_transcript_INTRON.bed

#wd=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/mergedBW

bed12s=( $prom ${utr5} ${cds} ${utr3} ${tes} )
bodylens=( 247 46 245 305 247 )
regions=( "prom" "utr5" "cds" "utr3" "tes" )

	cd ${wd}
     	echo "Processing ${f}"
     	basef=$(basename "${f}" .bw)


	for ((j=0; j<${#bed12s[*]}; j++)); do
		bed12=${bed12s[j]}
		bodylen=${bodylens[j]}
		region=${regions[j]}

		echo "Processing region ${region}"
		
		
			computeMatrix scale-regions \
			-S ${f} \
			-R ${bed12} \
			-b 0 \
			-a 0 \
			--metagene \
			--binSize 1 \
			-p max/2 \
			--missingDataAsZero \
			--regionBodyLength ${bodylen} \
			-o ${basef}_over_${region}_clip_matrix.gz
		
	done
		
	 if [ 1 -eq 0 ]; then
		echo "procesing intron"
		computeMatrix scale-regions \
		-S ${f} \
		-R ${intron} \
		-b 0 \
		-a 0 \
		--binSize 1 \
		--missingDataAsZero \
		--regionBodyLength 500 \
		-p max/2 \
		-o ${basef}_over_intron_clip_matrix.gz
		
		plotProfile \
		-m ${basef}_over_intron_clip_matrix.gz \
		-o ${outd}/${basef}_over_intron_clip.png \
		--plotType se \
		--averageType mean \
                --startLabel 5-prime \
                --endLabel 3-prime \
		--colors blue red green cyan brown black \
		--regionsLabel intron \
		--samplesLabel ${basef} \
		--numPlotsPerRow 1 \
		--perGroup # grouping by regions

	fi
