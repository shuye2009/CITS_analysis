#!/bin/bash
# This pipeline is used after CITS jobs have been finished 

TARGETbed=$1 ## name of CLIP target  
TARGETdir=$2 ## /home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/CLIP_cits/CITS_*  ## the location of the cits file
INPUTbed=$3 ## name of CLIP target
INPUTdir=$4 ## /home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/CLIP_cits/CITS_*  ## the location of the cits file
scriptDir=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/scripts_program/CITS_CIMS
plotSwitch=$5  ## "1,1,1,1,1,1
OUTdir=$6
mkdir -p $OUTdir

switches=($(echo $plotSwitch | tr "," " "))

# plot metaplotR for individual target, segmented in UTR5, CDS, UTR3, TES+2k
if [ 1 -eq ${switches[0]} ]; then
	echo "1.1 computing distance matrices for metaPlotR"
#	if [ -f $TARGETdir/${TARGETbed}.dist.measures.txt ]; then
		processing_input_bed.sh ${TARGETbed} ${TARGETdir} 
#	fi
	echo "1.2 producing metaPlotR plot"
	Rscript ${scriptDir}/metaPlotR_CLIP.R ${TARGETbed}.dist.measures.txt $TARGETdir $OUTdir
fi

# plot metagene for individual target, segmented in UTR5, CDS, UTR3, TES+2k
if [ 1 -eq ${switches[1]} ]; then
	echo "2.1 computing metagene profile"
	
#	if [ -f $TARGETdir/${TARGETbed}.bw ]; then
		cd $TARGETdir
		sort -k1,1 -k2,2n ${TARGETbed} | bedtools genomecov -i stdin -bg -g $HOME/genomic_feature/human.hg19.genome > ${TARGETbed}.bg
		bedGraphToBigWig ${TARGETbed}.bg $HOME/genomic_feature/human.hg19.genome ${TARGETbed}.bw
		${scriptDir}/deepTools_computeMatrix_metagene_parts.sh ${TARGETbed}.bw $TARGETdir $TARGETdir
#	fi

        echo "2.2 plot metagene profile"
        Rscript ${scriptDir}/metagene_plot_of_CITSpeak.R $TARGETbed $TARGETdir $OUTdir
fi

## only use this section for run_plot_CIMS_CITS_combined_proteinCoding_jobs.sh, plot from TSS to TES
if [ 1 -eq ${switches[2]} ]; then
        echo "3 plotting protein-coding gene profile"
        cd $TARGETdir
#	if [ -f ${TARGETbed}.bw ]; then
        	sort -k1,1 -k2,2n ${TARGETbed} | bedtools genomecov -i stdin -bg -g $HOME/genomic_feature/human.hg19.genome > ${TARGETbed}.bg
        	bedGraphToBigWig ${TARGETbed}.bg $HOME/genomic_feature/human.hg19.genome ${TARGETbed}.bw
#	fi
        ${scriptDir}/deepTools_computeMatrix_metagene_bamCompare_bc.sh ${TARGETbed}.bw $TARGETdir $OUTdir
fi

## only use section for run_plot_CIMS_CITS_combined_siSP1_jobs.sh, plot PAU up, down, noChange genes
if [ 1 -eq ${switches[3]} ]; then
        echo "4 plotting siSP1 PAU gene profile"
        cd $TARGETdir
 #       if [ -f ${TARGETbed}.bw ]; then
                sort -k1,1 -k2,2n ${TARGETbed} | bedtools genomecov -i stdin -bg -g $HOME/genomic_feature/human.hg19.genome > ${TARGETbed}.bg
                bedGraphToBigWig ${TARGETbed}.bg $HOME/genomic_feature/human.hg19.genome ${TARGETbed}.bw
 #      fi
        ${scriptDir}/deepTools_computeMatrix_siSP1_peak_profile.sh siSP1 false ${TARGETbed}.bw $TARGETdir $OUTdir
fi


############ TARGET Input overlay plot ####################

if [ 1 -eq ${switches[4]} ]; then
        echo "5. plot target-input metagene profile with metaPlotR_CLIP_Input.R"
	if [ ! -f ${INPUTdir}/${INPUTbed}.dist.measures.txt ]; then
		processing_input_bed.sh ${INPUTbed} ${INPUTdir}
	fi
#	if [ -f ${TARGETdir}/${TARGETbed}.dist.measures.txt ]; then
                processing_input_bed.sh ${TARGETbed} ${TARGETdir}
#       fi
        Rscript ${scriptDir}/metaPlotR_CLIP_Input.R $TARGETbed,$INPUTbed $TARGETdir,$INPUTdir $OUTdir FALSE
fi

if [ 1 -eq ${switches[5]} ]; then
        echo "6. plot target-input metagene profile with metagene_plot_of_CITSpeak.R"

#	if [ -f $TARGETdir/${TARGETbed}.bw ]; then
                cd $TARGETdir
                sort -k1,1 -k2,2n ${TARGETbed} | bedtools genomecov -i stdin -bg -g $HOME/genomic_feature/human.hg19.genome > ${TARGETbed}.bg
                bedGraphToBigWig ${TARGETbed}.bg $HOME/genomic_feature/human.hg19.genome ${TARGETbed}.bw
		${scriptDir}/deepTools_computeMatrix_metagene_parts.sh ${TARGETbed}.bw $TARGETdir $TARGETdir
#        fi

#	if [ -f $INPUTdir/${INPUTbed}.bw ]; then
                cd $INPUTdir
                sort -k1,1 -k2,2n ${INPUTbed} | bedtools genomecov -i stdin -bg -g $HOME/genomic_feature/human.hg19.genome > ${INPUTbed}.bg
                bedGraphToBigWig ${INPUTbed}.bg $HOME/genomic_feature/human.hg19.genome ${INPUTbed}.bw
		${scriptDir}/deepTools_computeMatrix_metagene_parts.sh ${INPUTbed}.bw $INPUTdir $INPUTdir
#       fi

        Rscript ${scriptDir}/metagene_plot_of_CITSpeak.R $TARGETbed,$INPUTbed $TARGETdir,$INPUTdir $OUTdir
fi


echo "Finish all"
