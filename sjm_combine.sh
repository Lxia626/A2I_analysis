#SAMPLE_ID=(SRR3192400 SRR3192401 SRR3192403 SRR3192406) ### need to modify according to each dataset
#SJM_list=()
#a=0
#for SAMPLE in ${SAMPLE_ID[@]}; do
#    SJM_list[a]=${SAMPLE}_lxia_star/${SAMPLE}_1_star.sjm
#    SAMPLE_list[a]=${SAMPLE}_1
#    a=`expr ${a} + 1`
#    SJM_list[a]=${SAMPLE}_lxia_star/${SAMPLE}_2_star.sjm
#    SAMPLE_list[a]=${SAMPLE}_2
#    a=`expr ${a} + 1`
#done
## echo ${SJM_list[@]}
## echo ${SAMPLE_list[@]}
#if [[ ! -d `pwd`/logs_pipeline01 ]]; then mkdir `pwd`/logs_pipeline01; fi

unique_id=$1 # need to specify "site01" or "site02" in ${unique_id}
if [[ -f ${unique_id}.order ]]; then rm ${unique_id}.order; fi
if [[ -f ${unique_id}.cmd ]]; then rm ${unique_id}.cmd; fi
if [[ -f ${unique_id}.sjm ]]; then rm ${unique_id}.sjm; fi

## step 1:
x=0
for SAMPLE in ${SAMPLE_ID[@]}; do
    sample_name=${SAMPLE_ID[x]}
    grep 'order' ${SAMPLE}/${sample_name}_1_star.sjm >> ${unique_id}.order
    grep 'order' ${SAMPLE}/${sample_name}_2_star.sjm >> ${unique_id}.order
    x=`expr ${x} + 1`
    echo -e "order ${sample_name}_1_GATKCallEdits before ${SAMPLE}_CHANGE" >> ${unique_id}.order # ${SAMPLE_list[x]}_First
    echo -e "order ${sample_name}_2_GATKCallEdits before ${SAMPLE}_CHANGE" >> ${unique_id}.order
    if [[ ${x} -lt ${#SAMPLE_ID[@]} ]]; then
        echo -e "order ${SAMPLE}_CHANGE before ${SAMPLE_ID[x]}_1_First" >> ${unique_id}.order
        echo -e "order ${SAMPLE}_CHANGE before ${SAMPLE_ID[x]}_2_First" >> ${unique_id}.order
    fi
    grep -vE 'order|log_dir' ${SAMPLE}/${sample_name}_1_star.sjm >> ${unique_id}.cmd
    grep -vE 'order|log_dir' ${SAMPLE}/${sample_name}_2_star.sjm >> ${unique_id}.cmd    
    echo -e "job_begin\n name ${SAMPLE}_CHANGE\n memory 3G\n directory `pwd`\n cmd_begin" >> ${unique_id}.cmd
    echo -e "/usr/bin/date;" >> ${unique_id}.cmd
    echo -e "cmd_end\njob_end\n" >> ${unique_id}.cmd
done

echo "log_dir `pwd`/logs_pipeline01" >> ${unique_id}.order
cat ${unique_id}.cmd ${unique_id}.order > ${unique_id}.sjm
# sjm ${unique_id}.sjm

### step 2:
x=0
Alu=$2
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
# sjm ${unique_id}.sjm
