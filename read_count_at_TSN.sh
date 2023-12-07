#!/bin/bash

## count number of clip reads at the transicription start nucleotide

txfile=$HOME/genomic_feature/gencode.v19.transcript_protein_coding.bed
utrfile=$HOME/genomic_feature/gencode.v19.protein_coding_transcript_UTR5.bed12

size=20  # size of extension
slopped=$HOME/genomic_feature/gencode.v19.transcript_protein_coding_TSN_${size}b.bed

## get the TSN, i.e. the first nucleotide at transcription start site, then extend $size bases on both sides
awk '{if($6 == "+") {print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$6} else {print $1"\t"$3-1"\t"$3"\t"$4"\t"$5"\t"$6}}' $utrfile | bedtools slop -b $size -g $HOME/genomic_feature/human_hg19.genome -i stdin >  $slopped

## get strand-aware coverage of intervals in a by intervals in b, output tag count if the coverage is greater than 0 
wd=/home/greenblattlab/shuyepu/Nabeel/CLIP2021/CLIP_cits/CITS_m6A
cd $wd

for tag in m6A_*.thUni.tag.bed; do
	bedtools coverage -s -c -a $slopped -b $tag | awk '{if($7 > 0 ) {print $0}}' > $tag.TSN.coverage
	bedtools coverage -s -c -a $utrfile -b $tag | awk '{if($7 > 0 ) {print $0}}' > $tag.5utr.coverage
done


