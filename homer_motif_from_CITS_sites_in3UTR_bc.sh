#!/bin/bash
# define motif use meme from piranha peak file
factor=$1
wd=$2
outd=$3 # /home/greenblattlab/shuyepu/Nabeel/clip_analysis/results_output/homer_analysis
mkdir -p ${outd}

Inputdir=$4 # /home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/Input_cits
cd ${wd}
f=$5

genome=~/hg19/human.hg19_sorted.genome
fasta=~/hg19/hg19.fa

#for f in final_combined_${factor}_crosslink_sites.bed; do 
#for f in combined_crosslink_${factor}.bed.recurrent_in_input.bed; do	

##f=${factor}_CITS_crosslinkSite_in_3UTR.bed	
echo -e "\n\nProcessing ${f}"
basef=$(basename "${f}" .bed)

extbed=${basef}_shortExtended.bed
	
echo "Extending bed ..."
bedtools slop -s -l 5 -r 5 -g ${genome} -i ${f} > ${extbed}
echo "Getting DNA sequence ..."
bedtools getfasta -fi ${fasta} -bed ${extbed} -s -fo ${extbed}.fa

motifd=${outd}/${basef}_ext5_homer_motif4_6_in3UTR
mkdir -p ${motifd}
	
echo "Running HOMER..."

#if test "$factor" = "Input"; then
	### use scrambled input sequence as background  
	scrambleFasta.pl ${extbed}.fa > ${factor}_3UTR_background.fa
	findMotifs.pl ${wd}/${extbed}.fa fasta ${motifd} -fasta ${factor}_3UTR_background.fa -rna -S 10 -mis 1 -p 5 -len 4,5,6 -basic -noknown
#else
	### use real input peak sequence as background
#	findMotifs.pl ${wd}/${extbed}.fa fasta ${motifd} -fasta ${Inputdir}/Input_CITS_crosslinkSite_in_3UTR_shortExtended.bed.fa -rna -S 10 -mis 1 -p 5 -len 4,5,6 -basic -noknown
#fi

#### use genome as background
#submitjob 100 -m 5 findMotifsGenome.pl $f hg19 ${outd}/${factor}_CITS_rnaMotif4_6 -size 15 -rna -S 10 -mis 1 -p 4 -len 4,5,6 -basic -noknown
#findMotifsGenome.pl $f hg19 ${outd}/${factor}_CITS_rnaMotif4_6 -size 15 -rna -S 10 -mis 1 -p 4 -len 4,5,6 -basic -noknown
	
echo "Running MEME..."
memed=${outd}/${basef}_ext5_meme_motif
#meme-chip ${extbed}.fa -dna -meme-nmotifs 5 -meme-mod anr -meme-p 8 -dreme-e 0.05 -meme-minw 6 -meme-maxw 5 -meme-minsites 100 -meme-maxsites 50 -oc ${memed}	
#meme ${extbed}.fa -dna -nmotifs 5 -mod zoops -p 8 -evt 0.05 -minw 6 -maxw 5 -minsites 100 -oc ${memed} -maxsize 10000000
#meme ${extbed}.fa -dna -nmotifs 5 -mod oops -p 4 -minw 4 -maxw 8 -oc ${memed} -maxsize 10000000
#dreme -p ${extbed}.fa -dna -m 3 -v 5 -e 0.05 -oc ${memed} -png -g 1000




