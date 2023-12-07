#! /bin/bash

cmd=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/scripts_program/CITS_CIMS/CITS_analysis.sh

# test if the CITS_analysis.sh run seccessfully, if not, move the bam fils to Zbam and run CITS analysis again on them
# for d in $(ls -d CITS_*); do cd $d; fa=${d//CITS_/}; if [[ ! -s "${fa}.tag.CITS.bed" ]]; then echo "$d $fa"; mv ../../CLIP_bam/$fa.bam* ../../CLIP_bam/Zbam/; fi; cd ..; done
# for d in $(ls -d CITS_*); do cd $d; fa=${d//CITS_/}; if [[ ! -s "${fa}.tag.CITS.bed" ]]; then echo "$d $fa"; fi; cd ..; done


## for CLIP2022 bams
if [ 1 -eq 0 ]; then
wd=/home/greenblattlab/shuyepu/Nabeel/CLIP2022/Input_bam
outDir=/home/greenblattlab/shuyepu/Nabeel/CLIP2022/CLIP_cits
mkdir -p $outDir

cd $wd

for f in *.bam; do
        echo "processing ${f}"
        basef=$(basename "${f}" .bam)
        outd=${outDir}/CITS_${basef}
        #if [[ ! -e $outd ]]; then
                echo "submitting ${f}"
                submitjob -w 200 -c 2 -m 20 $cmd $f $wd $outd
        #fi
done

fi


## for all CLIP bams
if [ 1 -eq 1 ]; then
wd=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/CLIP_bam
outDir=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/CLIP_cits
mkdir -p $outDir

cd $wd

for gene in KLF17 ZNF320 ZNF594 ZNF322 HuR RPS15A LSM6 ZNF23 InputPTBP1; do
        echo "processing ${gene}"
        for f in ${gene}.*.thUni.bam; do
	        echo "processing ${f}"
                basef=$(basename "${f}" .bam)
                outd=${outDir}/CITS_${basef}
        #if [[ ! -e $outd ]]; then   
       		echo "submitting ${f}"
	     	submitjob -w 200 -c 4 -m 20 $cmd $f $wd $outd
        #fi
        done
done

fi


## for CLIP2021
if [ 1 -eq 0 ]; then
wd=/home/greenblattlab/shuyepu/Nabeel/CLIP2021/CLIP_bam
outDir=/home/greenblattlab/shuyepu/Nabeel/CLIP2021/CLIP_cits
mkdir -p $outDir

cd $wd

for f in NXF1.B.thUni.bam YTHDF1_A.thUni.bam ZNF384.A.thUni.bam ADAR1_C.thUni.bam CHTOP.A.thUni.bam MATR3_A.thUni.bam Alyref_S2.thUni.bam; do
	if [[ ! $f =~ .*m6A.* ]]; then 
                echo "processing ${f}"
		basef=$(basename "${f}" .bam)
                outd=${outDir}/CITS_${basef}
                submitjob 100 -m 50 $cmd $f $wd $outd
        fi

done

fi

## for all Input
if [ 1 -eq 0 ]; then
wd=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/Input_bam
outDir=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/CLIP_cits
mkdir -p $outDir

cd $wd

for factor in $(ls *.bam | cut -d. -f1 | sort | uniq | tr "\n" " "); do
#for factor in SP1_S2; do
        for f in ${factor}*.bam; do
                echo "processing ${f}"
		basef=$(basename "${f}" .bam)
                outd=${outDir}/CITS_${basef}
                submitjob 200 -m 40 $cmd $f $wd $outd
        done

done

fi


## for SP1_2021
if [ 1 -eq 0 ]; then
wd=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/SP1_2021/CLIP_bam
outDir=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/SP1_2021/CLIP_cits
mkdir -p $outDir

cd $wd

#for factor in $(ls *.bam | cut -d. -f2 | sort | uniq | tr "\n" " "); do
for factor in SP1_S2; do
	for f in *.${factor}.*.bam; do
        	echo "processing ${f}"
        	basef=$(basename "${f}" .bam)
        	outd=${outDir}/CITS_${factor}
        	submitjob 200 -m 50 $cmd $f $wd $outd 0.05
	done

done 
fi

## for ac4c_clip
if [ 1 -eq 0 ]; then
wd=/home/greenblattlab/shuyepu/Nabeel/ac4c_clip
outDir=/home/greenblattlab/shuyepu/Nabeel/ac4c_clip/CITS_CIMS
mkdir -p $outDir

cd $wd

#for factor in $(ls *.bam | cut -d. -f2 | sort | uniq | tr "\n" " "); do
for factor in Ac4c; do
        for f in *_${factor}_*.bam; do
                echo "processing ${f}"
                basef=$(basename "${f}" .bam)
                outd=${outDir}/${basef}
                submitjob 200 -m 50 $cmd $f $wd $outd 0.05
        done

done
fi


