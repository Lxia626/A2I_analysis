#! /usr/bin/bash
 
STAR_INDEX=$1 # /public/home/Songlab/xial/Work/cell_line_dataset/HEK293T/GSE99249/STAR_index/
UNIQUE_ID=$2 # GSE99249_NF

# step 1
#ls -l | grep '^d' | awk '{print $NF}' | while read SAMPLE; do
#    cd ${SAMPLE}
#    sh ~/scripts/pipeline/A2I_editing_pipe/site_01_readsFiltrate_star.sh `pwd` ${STAR_INDEX} &
#    cd ../
#done
#
#mkdir logs_pipeline01
#unique_id="${UNIQUE_ID}_site01" # not filter
#if [[ -f ${unique_id}.order ]]; then rm ${unique_id}.order; fi
#if [[ -f ${unique_id}.cmd ]]; then rm ${unique_id}.cmd; fi
#if [[ -f ${unique_id}.sjm ]]; then rm ${unique_id}.sjm; fi
#SAMPLE_ID=($(ls -l | grep '^d' | awk '{print $NF}'))
#x=0
#for SAMPLE in ${SAMPLE_ID[@]}; do
#    sample_name=${SAMPLE_ID[x]}
#    grep 'order' ${SAMPLE}/${sample_name}_1_star.sjm >> ${unique_id}.order
#    grep 'order' ${SAMPLE}/${sample_name}_2_star.sjm >> ${unique_id}.order
#    x=`expr ${x} + 1`
#    echo -e "order ${sample_name}_1_GATKCallEdits before ${SAMPLE}_CHANGE" >> ${unique_id}.order # ${SAMPLE_list[x]}_First
#    echo -e "order ${sample_name}_2_GATKCallEdits before ${SAMPLE}_CHANGE" >> ${unique_id}.order
#    if [[ ${x} -lt ${#SAMPLE_ID[@]} ]]; then
#        echo -e "order ${SAMPLE}_CHANGE before ${SAMPLE_ID[x]}_1_First" >> ${unique_id}.order
#        echo -e "order ${SAMPLE}_CHANGE before ${SAMPLE_ID[x]}_2_First" >> ${unique_id}.order
#    fi
#    grep -vE 'order|log_dir' ${SAMPLE}/${sample_name}_1_star.sjm >> ${unique_id}.cmd
#    grep -vE 'order|log_dir' ${SAMPLE}/${sample_name}_2_star.sjm >> ${unique_id}.cmd
#    echo -e "job_begin\n name ${SAMPLE}_CHANGE\n memory 3G\n directory `pwd`\n cmd_begin" >> ${unique_id}.cmd
#    echo -e "/usr/bin/date;" >> ${unique_id}.cmd
#    echo -e "cmd_end\njob_end\n" >> ${unique_id}.cmd
#done
#echo "log_dir `pwd`/logs_pipeline01" >> ${unique_id}.order
#cat ${unique_id}.cmd ${unique_id}.order > ${unique_id}.sjm

# step 2
ls -l | grep '^d' | awk '{print $NF}' | while read SAMPLE; do
    cd ${SAMPLE}
    sh ~/scripts/pipeline/A2I_editing_pipe/site_02_callSite.sh `pwd` Both ${SAMPLE}_NF
    cd ../
done

mkdir logs_pipeline02
unique_id="${UNIQUE_ID}_site02"
SAMPLE_ID=($(ls -l | grep '^d' | awk '{print $NF}'))
Alu="Both"
for SAMPLE in ${SAMPLE_ID[@]}; do
    sample_name=${SAMPLE_ID[x]}
    grep 'order' ${SAMPLE}/${sample_name}.bam.sjm >> ${unique_id}.order
    x=`expr ${x} + 1`
    if [[ ${Alu} == "Alu" || ${Alu} == "nonAlu" ]]; then
        echo -e "order ${sample_name}_CAT before ${SAMPLE}_CHANGE" >> ${unique_id}.order
    else
        echo -e "order ${sample_name}_CAT_NoAlu before ${SAMPLE}_CHANGE" >> ${unique_id}.order
    fi
    if [[ ${x} -lt ${#SAMPLE_ID[@]} ]]; then
        echo -e "order ${SAMPLE}_CHANGE before ${SAMPLE_ID[x]}_MERGE" >> ${unique_id}.order
    fi
    grep -vE 'order|log_dir' ${SAMPLE}/${sample_name}.bam.sjm >> ${unique_id}.cmd
    echo -e "job_begin\n name ${SAMPLE}_CHANGE\n memory 3G\n directory `pwd`\n cmd_begin" >> ${unique_id}.cmd
    echo -e "/usr/bin/date;" >> ${unique_id}.cmd
    echo -e "cmd_end\njob_end\n" >> ${unique_id}.cmd
done
echo "log_dir `pwd`/logs_pipeline02" >> ${unique_id}.order
cat ${unique_id}.cmd ${unique_id}.order > ${unique_id}.sjm
