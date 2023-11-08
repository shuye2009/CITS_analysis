#!/bin/bash

basef=$1
outd=$2
pv=$3

# for Nabeel CLIP pv=0.05, for YTHD downlaoded CLIP pv=0.0001
# pv=0.0001

## for CIMS use tag count >= 10 and mutation count >= 5 and FDR <= $pv
if [[ -e ${outd}/${basef}.tag.sub.CIMS.txt ]]; then
	awk '{if($9<=$pv) {print $0}}' ${outd}/${basef}.tag.sub.CIMS.txt | sort -k 9,9n -k 8,8nr -k 7,7n | cut -f 1-6 > ${outd}/${basef}.tag.sub.CIMS.${pv}.bed
fi

if [[ -e ${outd}/${basef}.tag.del.CIMS.txt ]]; then
	awk '{if($9<=$pv) {print $0}}' ${outd}/${basef}.tag.del.CIMS.txt | sort -k 9,9n -k 8,8nr -k 7,7n | cut -f 1-6 > ${outd}/${basef}.tag.del.CIMS.${pv}.bed
fi

$HOME/ctk-1.1.4/thresholdCITS.pl ${outd}/${basef}.tag.CITS.bed ${outd}/${basef}.tag.CITS.${pv}.bed $pv

#cat ${outd}/${basef}.tag.del.CIMS.bed ${outd}/${basef}.tag.sub.CIMS.bed > ${outd}/${basef}.tag.CIMS.bed
cat ${outd}/${basef}.tag.del.CIMS.${pv}.bed ${outd}/${basef}.tag.sub.CIMS.${pv}.bed > ${outd}/${basef}.tag.CIMS.${pv}.bed
cat ${outd}/${basef}.tag.CITS.${pv}.bed ${outd}/${basef}.tag.CIMS.${pv}.bed > ${outd}/${basef}.crossLink_site.${pv}.bed
