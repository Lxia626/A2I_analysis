#! /usr/bin/bash

# STAR_INDEX:
reads_length=xxx # need to specify accordingly
mkdir STAR_INDEX
cd STAR_INDEX
STAR --runThreadN 10 --sjdbOverhang `expr ${reads_length} - 1` --runMode genomeGenerate --sjdbGTFfile /public/home/Songlab/songyl/database/annotation/human/gencode.v19.chr_patch_hapl_scaff.annotation.gtf --genomeFastaFiles /public/home/Songlab/songyl/database/genome/human/GencodeV19/STAR_index/with_GencodeV19_gtf/GRCh37.p13.genome.fa --genomeDir ./ > STAR_index.log 2>&1 &
cd ../
wait

# Basic options
FASTQ_DIR=xxx # need to specify accordingly
home_wd=`pwd`
STAR_INDEX=${home_wd}/STAR_INDEX # /need to specify accordingly, e.g. /public/home/Songlab/xial/Work/cell_line_dataset/HEK293T/GSE99249/STAR_index/
UNIQUE_ID=${home_wd##*/} # to specify if needed
SraFile=FASTQ/SraAccList.csv # to specify if needed
SAMPLEs=($(cat ${SraFile}))
Alu='Both' # to specify if needed {'Both', 'Alu', 'nonAlu'}

# construct directory structure
for sample in ${SAMPLEs[@]}; do
    mkdir ${sample}
    cd ${sample}
    for FASTQ in `ls ${FASTQ_DIR}/${sample}*`; do
        ln -s ${FASTQ} .
    done
    cd ../
done

mkdir logs_pipeline
# step 1
for SAMPLE in ${SAMPLEs[@]}; do
    cd ${SAMPLE}
    sh ~/scripts/pipeline/A2I_editing_pipe/A2I_pipe/site_01_readsFiltrate_star.sh `pwd` ${STAR_INDEX} &
    cd ../
done

wait
if [[ -f ${UNIQUE_ID}.order ]]; then rm ${UNIQUE_ID}.order; fi
if [[ -f ${UNIQUE_ID}.cmd ]];   then rm ${UNIQUE_ID}.cmd; fi
if [[ -f ${UNIQUE_ID}.sjm ]];   then rm ${UNIQUE_ID}.sjm; fi
x=0
for SAMPLE in ${SAMPLEs[@]}; do
    sample_name=${SAMPLEs[x]}
    grep 'order' ${SAMPLE}/${sample_name}_1_star.sjm >> ${UNIQUE_ID}.order
    grep 'order' ${SAMPLE}/${sample_name}_2_star.sjm >> ${UNIQUE_ID}.order
    touch ${SAMPLE}/starsplice_${SAMPLE}_1/${sample_name}_1_GATKRECAL.filtered.bam
    touch ${SAMPLE}/starsplice_${SAMPLE}_2/${sample_name}_2_GATKRECAL.filtered.bam
    x=`expr ${x} + 1`
    echo -e "order ${sample_name}_1_GATKCallEdits before ${SAMPLE}_CHANGE" >> ${UNIQUE_ID}.order # ${SAMPLE_list[x]}_First
    echo -e "order ${sample_name}_2_GATKCallEdits before ${SAMPLE}_CHANGE" >> ${UNIQUE_ID}.order
    if [[ ${x} -lt ${#SAMPLEs[@]} ]]; then
        echo -e "order ${SAMPLE}_CHANGE before ${SAMPLEs[x]}_1_First" >> ${UNIQUE_ID}.order
        echo -e "order ${SAMPLE}_CHANGE before ${SAMPLEs[x]}_2_First" >> ${UNIQUE_ID}.order
    fi
    grep -vE 'order|log_dir' ${SAMPLE}/${sample_name}_1_star.sjm >> ${UNIQUE_ID}.cmd
    grep -vE 'order|log_dir' ${SAMPLE}/${sample_name}_2_star.sjm >> ${UNIQUE_ID}.cmd
    echo -e "job_begin\n name ${SAMPLE}_CHANGE\n memory 3G\n directory `pwd`\n cmd_begin" >> ${UNIQUE_ID}.cmd
    echo -e "/usr/bin/date;" >> ${UNIQUE_ID}.cmd
    echo -e "cmd_end\njob_end\n" >> ${UNIQUE_ID}.cmd
done
echo -e "job_begin\n name site01_CHANGE\n memory 3G\n directory `pwd`\n cmd_begin" >> ${UNIQUE_ID}.cmd
echo -e "/usr/bin/date;" >> ${UNIQUE_ID}.cmd
echo -e "cmd_end\njob_end\n" >> ${UNIQUE_ID}.cmd

# step 2
for SAMPLE in ${SAMPLEs[@]}; do
    cd ${SAMPLE}
    sh ~/scripts/pipeline/A2I_editing_pipe/A2I_pipe/site_02_callSite.sh `pwd` Both ${SAMPLE}
    cd ../
done

x=0
for SAMPLE in ${SAMPLEs[@]}; do
    sample_name=${SAMPLEs[x]}
    echo -e "order site01_CHANGE before ${sample_name}_MERGE" >> ${UNIQUE_ID}.order
    grep 'order' ${SAMPLE}/${sample_name}.bam.sjm >> ${UNIQUE_ID}.order
    x=`expr ${x} + 1`
    if [[ ${Alu} == "Alu" || ${Alu} == "nonAlu" ]]; then
        echo -e "order ${sample_name}_CAT before ${SAMPLE}_CHANGE" >> ${UNIQUE_ID}.order
    else
        echo -e "order ${sample_name}_CAT_NoAlu before ${SAMPLE}_CHANGE" >> ${UNIQUE_ID}.order
    fi
    if [[ ${x} -lt ${#SAMPLEs[@]} ]]; then
        echo -e "order ${SAMPLE}_CHANGE before ${SAMPLEs[x]}_MERGE" >> ${UNIQUE_ID}.order
    fi
    grep -vE 'order|log_dir' ${SAMPLE}/${sample_name}.bam.sjm >> ${UNIQUE_ID}.cmd
    echo -e "job_begin\n name ${SAMPLE}_CHANGE\n memory 3G\n directory `pwd`\n cmd_begin" >> ${UNIQUE_ID}.cmd
    echo -e "/usr/bin/date;" >> ${UNIQUE_ID}.cmd
    echo -e "cmd_end\njob_end\n" >> ${UNIQUE_ID}.cmd
done
echo "log_dir `pwd`/logs_pipeline" >> ${UNIQUE_ID}.order
cat ${UNIQUE_ID}.cmd ${UNIQUE_ID}.order > ${UNIQUE_ID}.sjm
