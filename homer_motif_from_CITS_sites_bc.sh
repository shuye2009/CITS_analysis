#!/bin/bash
# define motif use homer from CITS peak file
factor=$1
wd=$2
outd=$3  # /home/greenblattlab/shuyepu/Nabeel/clip_analysis/results_output/homer_analysis
Inputdir=$4  # /home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/Input_cits
f=$5  # bed file from which motif is to be found

cd ${wd}
mkdir -p ${outd}

genome=~/hg19/human.hg19_sorted.genome
fasta=~/hg19/hg19.fa

#for f in final_combined_${factor}_crosslink_sites.bed; do 
#for f in combined_crosslink_${factor}.bed.recurrent_in_input.bed; do	

##f=${factor}_over_Input_DESeq2_significant_crossLink_site.bed	
echo -e "\n\nProcessing ${f}"
basef=$(basename "${f}" .bed)

extbed=${basef}_shortExtended.bed
	
echo "Extending bed ..."
bedtools slop -s -l 5 -r 5 -g ${genome} -i ${f} > ${extbed}
echo "Getting DNA sequence ..."
bedtools getfasta -fi ${fasta} -bed ${extbed} -s -fo ${extbed}.fa

motifd=${outd}/${basef}_ext5_homer_motif4_6_allSites
mkdir -p ${motifd}
	
echo "Running HOMER..."

#if test "$factor" = "Input"; then
	### use scrambled input sequence as background  
	scrambleFasta.pl ${extbed}.fa > ${factor}_background.fa
	findMotifs.pl ${wd}/${extbed}.fa fasta ${motifd} -fasta ${factor}_background.fa -rna -S 10 -mis 1 -p 5 -len 4,5,6 -basic -noknown
#else
	### use real input peak sequence as background
#	findMotifs.pl ${wd}/${extbed}.fa fasta ${motifd} -fasta ${Inputdir}/combined_crosslink_Input.uniq_shortExtended.bed.fa -rna -S 10 -mis 1 -p 5 -len 4,5,6 -basic -noknown
#fi

#### use genome as background
#submitjob 100 -m 5 findMotifsGenome.pl $f hg19 ${outd}/${factor}_CITS_rnaMotif4_6 -size 15 -rna -S 10 -mis 1 -p 4 -len 4,5,6 -basic -noknown
#findMotifsGenome.pl $f hg19 ${outd}/${factor}_CITS_rnaMotif4_6 -size 15 -rna -S 10 -mis 1 -p 4 -len 4,5,6 -basic -noknown
	
echo "Running MEME..."
memed=${outd}/${basef}_ext5_meme_motif
#meme ${extbed}.fa -dna -nmotifs 5 -mod oops -p 4 -minw 4 -maxw 8 -oc ${memed} -maxsize 10000000
#dreme -p ${extbed}.fa -dna -v 5 -norc -e 1e-100 -oc ${memed} -png -g 1000	




