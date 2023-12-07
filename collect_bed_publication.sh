#!/bin/bash

## for SP1 paper

cpsf="CPSF2.A CPSF2.B NUDT21.C NUDT21.D NUDT21.E CPSF7.A CPSF7.B CPSF160.A CPSF160.B FIP1.A FIP1.B"
cpsfdir="/home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/CLIP_cits" #CITS_*.thUni
sp1="SP1.A SP1.B"
sp1dir="/home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/SP1_analysis_AB/CLIP_cits"
#cpsf5="NUDT21.A NUDT21.B NUDT21.C NUDT21.D NUDT21.E"
#cpsf5dir=""
input="Input.CPA.A Input.NUDT21.A Input.SP1.A"
inputdir="/home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/Input_cits"

outdir="/home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/SP1_paper"
mkdir -p $outdir

for f in $cpsf; do
 filedir=${cpsfdir}/CITS_${f}.thUni
 outfile=${f}.thUni.crossLink_site.1e-10.bed
 cp $filedir/$outfile $outdir/$outfile
done

for f in $sp1; do
 filedir=${sp1dir}/CITS_${f}.thUni
 outfile=${f}.thUni.crossLink_site.1e-10.bed
 cp $filedir/$outfile $outdir/$outfile
done

for f in $input; do
 filedir=${inputdir}/CITS_${f}.thUni
 outfile=${f}.thUni.crossLink_site.1e-10.bed
 cp $filedir/$outfile $outdir/$outfile
done

## END


