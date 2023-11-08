#!/bin/bash

## find motif for 10 nt sequence around crosslink sites in the peak bed file
## use 10 nt sequence arround sites in second bed as backgound

readonly refsize=$HOME/GRCH37_gencode/GRCh37.p13.genome.fa.size
readonly reffasta=$HOME/GRCH37_gencode/GRCh37.p13.genome.fa

peakbed=$1
slop=$2
backbed=$3
alphabet=$4

if [[ -z $slop ]]; then slop=0; fi
if [[ -z $alphabet ]]; then alphabet="-rna"; fi

basef=$(basename $peakbed .bed)
echo $peakbed
echo $slop

if [[ ! -z $backbed ]]; then
	baseb=$(basename $backbed .bed)
	echo $backbed
	echo $baseb

	dremed=dreme_$basef
	mkdir -p $dremed
	bedtools slop -i $peakbed -g $refsize -b $slop | bedtools getfasta -bed stdin -s -name -fi $reffasta -fo ${basef}.fa
	bedtools slop -i $backbed -g $refsize -b $slop | bedtools getfasta -bed stdin -s -name -fi $reffasta -fo ${baseb}.fa

	dreme -p ${basef}.fa -n ${baseb}.fa -norc -k 6 -m 3 $alphabet -oc $dremed


else
	dremed=dreme_${basef}_noInput
	mkdir -p $dremed
	fl=$((slop*2+1))
	bedtools slop -i $peakbed -g $refsize -b $slop | bedtools getfasta -bed stdin -s -name -fi $reffasta -fo ${basef}.fa
	bedtools slop -i $peakbed -g $refsize -b $slop | bedtools flank -i stdin -g $refsize -b $fl | bedtools getfasta -bed stdin -s -name -fi $reffasta -fo ${basef}_bg.fa
	#meme ${basef}.fa -p 1 -neg ${basef}_bg.fa -nmotifs 3 -objfun de -minw 5 -maxw 7 -evt 0.01 -rna -mod zoops -oc ${dremed}
	dreme -p ${basef}.fa -n ${basef}_bg.fa -norc -mink 6 -maxk 6 -m 5 -e 0.01 $alphabet -oc ${dremed}
	#streme --p ${basef}.fa --n ${basef}_bg.fa --objfun de --nmotifs 5 --minw 5 --maxw 7 --pvt 0.01 --rna --oc ${dremed}
fi
## END
