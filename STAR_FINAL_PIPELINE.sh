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

# construct directory structure
for sample in ${SAMPLEs[@]}; do
    mkdir ${sample}
    cd ${sample}
    for FASTQ in `ls ${FASTQ_DIR}/${sample}*`; do
        ln -s ${FASTQ} .
    done
    cd ../
done

# step 1
for SAMPLE in ${SAMPLEs[@]}; do
    cd ${SAMPLE}
    sh ~/scripts/pipeline/A2I_editing_pipe/A2I_pipe/site_01_readsFiltrate_star.sh `pwd` ${STAR_INDEX} &
    cd ../
done

wait
mkdir logs_pipeline01
unique_id="${UNIQUE_ID}_site01" # not filter
if [[ -f ${unique_id}.order ]]; then rm ${unique_id}.order; fi
if [[ -f ${unique_id}.cmd ]]; then rm ${unique_id}.cmd; fi
if [[ -f ${unique_id}.sjm ]]; then rm ${unique_id}.sjm; fi
x=0
for SAMPLE in ${SAMPLEs[@]}; do
    sample_name=${SAMPLEs[x]}
    grep 'order' ${SAMPLE}/${sample_name}_1_star.sjm >> ${unique_id}.order
    grep 'order' ${SAMPLE}/${sample_name}_2_star.sjm >> ${unique_id}.order
    touch ${SAMPLE}/starsplice_${SAMPLE}_1/${sample_name}_1_GATKRECAL.filtered.bam
    touch ${SAMPLE}/starsplice_${SAMPLE}_2/${sample_name}_2_GATKRECAL.filtered.bam
    x=`expr ${x} + 1`
    echo -e "order ${sample_name}_1_GATKCallEdits before ${SAMPLE}_CHANGE" >> ${unique_id}.order # ${SAMPLE_list[x]}_First
    echo -e "order ${sample_name}_2_GATKCallEdits before ${SAMPLE}_CHANGE" >> ${unique_id}.order
    if [[ ${x} -lt ${#SAMPLEs[@]} ]]; then
        echo -e "order ${SAMPLE}_CHANGE before ${SAMPLEs[x]}_1_First" >> ${unique_id}.order
        echo -e "order ${SAMPLE}_CHANGE before ${SAMPLEs[x]}_2_First" >> ${unique_id}.order
    fi
    grep -vE 'order|log_dir' ${SAMPLE}/${sample_name}_1_star.sjm >> ${unique_id}.cmd
    grep -vE 'order|log_dir' ${SAMPLE}/${sample_name}_2_star.sjm >> ${unique_id}.cmd
    echo -e "job_begin\n name ${SAMPLE}_CHANGE\n memory 3G\n directory `pwd`\n cmd_begin" >> ${unique_id}.cmd
    echo -e "/usr/bin/date;" >> ${unique_id}.cmd
    echo -e "cmd_end\njob_end\n" >> ${unique_id}.cmd
done
#echo -e "job_begin\n name site01_CHANGE\n memory 3G\n directory `pwd`\n cmd_begin" >> ${unique_id}.cmd
#echo -e "/usr/bin/date;" >> ${unique_id}.cmd
#echo -e "cmd_end\njob_end\n" >> ${unique_id}.cmd
echo "log_dir `pwd`/logs_pipeline01" >> ${unique_id}.order
cat ${unique_id}.cmd ${unique_id}.order > ${unique_id}.sjm

# step 2
for SAMPLE in ${SAMPLEs[@]}; do
    cd ${SAMPLE}
    sh ~/scripts/pipeline/A2I_editing_pipe/A2I_pipe/site_02_callSite.sh `pwd` Both ${SAMPLE}
    cd ../
done

mkdir logs_pipeline02
unique_id="${UNIQUE_ID}_site02"
Alu="Both"
x=0
for SAMPLE in ${SAMPLEs[@]}; do
    sample_name=${SAMPLEs[x]}
    #echo -e "order site01_CHANGE before ${sample_name}_MERGE" >> ${unique_id}.order
    grep 'order' ${SAMPLE}/${sample_name}.bam.sjm >> ${unique_id}.order
    x=`expr ${x} + 1`
    if [[ ${Alu} == "Alu" || ${Alu} == "nonAlu" ]]; then
        echo -e "order ${sample_name}_CAT before ${SAMPLE}_CHANGE" >> ${unique_id}.order
    else
        echo -e "order ${sample_name}_CAT_NoAlu before ${SAMPLE}_CHANGE" >> ${unique_id}.order
    fi
    if [[ ${x} -lt ${#SAMPLEs[@]} ]]; then
        echo -e "order ${SAMPLE}_CHANGE before ${SAMPLEs[x]}_MERGE" >> ${unique_id}.order
    fi
    grep -vE 'order|log_dir' ${SAMPLE}/${sample_name}.bam.sjm >> ${unique_id}.cmd
    echo -e "job_begin\n name ${SAMPLE}_CHANGE\n memory 3G\n directory `pwd`\n cmd_begin" >> ${unique_id}.cmd
    echo -e "/usr/bin/date;" >> ${unique_id}.cmd
    echo -e "cmd_end\njob_end\n" >> ${unique_id}.cmd
done
echo "log_dir `pwd`/logs_pipeline02" >> ${unique_id}.order
cat ${unique_id}.cmd ${unique_id}.order > ${unique_id}.sjm
