#!/bin/bash

bamf=$1
inputd=$2
outd=$3
pv=0.05 ## default pvalue, can be thresholded later on

ctkdir=$HOME/ctk-1.1.4
if [ -d $outd ]; then rm -r $outd; fi
mkdir -p $outd

basef=$(basename "${bamf}" .bam)
if [[ ! -e $inputd/${basef}.sam ]]; then
	echo "bam to sam conversion"
	/usr/bin/samtools view --threads 2 $inputd/$bamf > $inputd/${basef}.sam
fi

perl $ctkdir/parseAlignment.pl -v --map-qual 40 --min-len 18 --mutation-file ${outd}/${basef}.mutation.txt $inputd/${basef}.sam ${outd}/${basef}.tag.bed
rm ${basef}.sam

python $ctkdir/joinWrapper.py ${outd}/${basef}.mutation.txt ${outd}/${basef}.tag.bed 4 4 N ${outd}/${basef}.tag.mutation.txt

perl $ctkdir/getMutationType.pl -t del ${outd}/${basef}.tag.mutation.txt ${outd}/${basef}.tag.del.bed
perl $ctkdir/getMutationType.pl -t ins ${outd}/${basef}.tag.mutation.txt ${outd}/${basef}.tag.ins.bed
perl $ctkdir/getMutationType.pl -t sub --summary ${outd}/${basef}.mutation_stats.txt ${outd}/${basef}.tag.mutation.txt ${outd}/${basef}.tag.sub.bed

mkdir /dev/shm/${basef}
perl $ctkdir/CIMS.pl -big -n 5 -p -outp ${outd}/${basef}.tag.sub.pos.stat.txt -v -FDR $pv -c /dev/shm/${basef}/cache_sub ${outd}/${basef}.tag.bed ${outd}/${basef}.tag.sub.bed ${outd}/${basef}.tag.sub.CIMS.txt
#awk '{if($9<=$pv) {print $0}}' ${outd}/${basef}.tag.sub.CIMS.txt | sort -k 9,9n -k 8,8nr -k 7,7n > ${outd}/${basef}.tag.sub.CIMS.${pv}.txt # significant at 3rd decimal position
cut -f 1-6 ${outd}/${basef}.tag.sub.CIMS.txt > ${outd}/${basef}.tag.sub.CIMS.bed
#awk '{print $1"\t"$2-10"\t"$3+10"\t"$4"\t"$5"\t"$6}' ${outd}/${basef}.tag.sub.CIMS.${pv}.bed > ${outd}/${basef}.tag.sub.CIMS.${pv}.21nt.bed
#+/-10 around CIMS

perl $ctkdir/CIMS.pl -big -n 5 -p -outp ${outd}/${basef}.tag.del.pos.stat.txt -v -FDR $pv -c /dev/shm/${basef}/cache_del ${outd}/${basef}.tag.bed ${outd}/${basef}.tag.del.bed ${outd}/${basef}.tag.del.CIMS.txt
#awk '{if($9<=$pv) {print $0}}' ${outd}/${basef}.tag.del.CIMS.txt | sort -k 9,9n -k 8,8nr -k 7,7n > ${outd}/${basef}.tag.del.CIMS.${pv}.txt # significant at 3rd decimal position
cut -f 1-6 ${outd}/${basef}.tag.del.CIMS.txt > ${outd}/${basef}.tag.del.CIMS.bed
#awk '{print $1"\t"$2-10"\t"$3+10"\t"$4"\t"$5"\t"$6}' ${outd}/${basef}.tag.del.CIMS.${pv}.bed > ${outd}/${basef}.tag.del.CIMS.${pv}.21nt.bed
#+/-10 around CIMS

perl $ctkdir/CITS.pl -big -p ${pv} --gap 25 -v -c /dev/shm/${basef}/cache_cits ${outd}/${basef}.tag.bed ${outd}/${basef}.tag.del.bed ${outd}/${basef}.tag.CITS.bed
#awk '{if($3-$2==1) {print $0}}' ${outd}/${basef}.tag.CITS.${pv}.bed > ${outd}/${basef}.tag.CITS.${pv}.singleton.bed
rm -r /dev/shm/${basef}

cat ${outd}/${basef}.tag.CITS.bed ${outd}/${basef}.tag.del.CIMS.bed ${outd}/${basef}.tag.sub.CIMS.bed > ${outd}/${basef}.crossLink_site.bed
#awk '{print $1"\t"$2-10"\t"$3+10"\t"$4"\t"$5"\t"$6}' ${outd}/${basef}.tag.CITS.${pv}.singleton.bed > ${outd}/${basef}.tag.CITS.${pv}.singleton.21nt.bed
