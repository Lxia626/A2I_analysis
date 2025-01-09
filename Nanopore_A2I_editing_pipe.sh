#! /usr/bin/env bash

# Setup getopt
help_info="Usage : \
\n\tsh $0 [options] \
\nOptions: \
\n\t-h, --help: to display help information \
\n\t-i, --input\tINPATH: Path where FASTQ files are stored. Note, the absolute path is recommended to be supplied \
\n\t-o, --output\tOUTPATH_ROOT: Path where to run the pipeline. Note, the absolute path is recommended to be supplied \
\n\t-u, --uniq\tUNIQ_ID: Unique ID for marking current run \
\n\t-L, --length\tLENGTH: Minimum length of reads [default: 0] \
\n\t-Q, --qual\tQUAL: Minimum base quality [default: 7] \
\n\t-d, --depth\tDEPTH: Minimum sequencing depth [default: 2] \
\n\t--suffix_ori\tSUFFIX_ORI: SUFFIX of raw FASTQ file [default: fastq] \
\n\t--multiPro\tMultiprocesses: Number of processing proceeding meanwhile [default: 30] \
\n\t--delete\tDELETE_ORIGINAL: Whether to delete the existed results in each output directory. 1: True, 0: False [default: 0] \
\n\t--level\tLEVEL: Whether to only calculate the editing level. 1:True, 0:False [default: 0]
\n\t
\nNote:
\n\tWhen running SJM, command 'sjm --max_running 15 \${sjmname}' is suggested"
ARGS=`getopt --options hi:o:u:L:Q:d: --long help,input:,output:,uniq:,length:,qual:,depth:,suffix_ori:,multiPro:,delete,level -n "$0" -- "$@"`
if [ $? != 0 ]; then 
    echo "Terminating..."
    exit
fi

echo ARGS=[$ARGS]
eval set -- "${ARGS}"
while true; do
    case "$1" in 
        -h|--help) echo -e ${help_info} >&2; exit 1; ;;
        -i|--input) INPATH=$2; shift; ;; # not need the absolute path
        -o|--output) OUTPATH_ROOT=$2; shift; ;; # not need the absolute path
        -u|--uniq) UNIQ_ID=$2; shift; ;; 
        -L|--length) LENGTH=$2; shift; ;; 
        -Q|--qual) QUAL=$2; shift; ;;
        -d|--depth) DEPTH=$2; shift; ;;
        --suffix_ori) SUFFIX_ORI=$2; shift; ;;
        --multiPro) MULTIPRO_NUM=$2; shift; ;;
        --delete) DELETE_ORIGINAL='1' ;; #{'0':FALSE, '1': TRUE}
        --level) LEVEL='1' ;; #{'0': FALSE, '1':TRUE} 
        :) echo "ERROR: -$OPTARG expects an corresponding argument" >&2; ;;
        --) shift; break; ;;
        *) echo "Internal error!"; exit 1; ;;
    esac
shift
done

if [ ! ${LENGTH} ]; then LENGTH=0; fi
if [ ! ${QUAL} ]; then QUAL=7; fi
if [ ! ${DEPTH} ]; then DEPTH=2; fi
if [ ! ${MULTIPRO_NUM} ]; then MULTIPRO_NUM=30; fi
if [ ! ${DELETE_ORIGINAL} ]; then DELETE_ORIGINAL='0'; fi
if [ ! ${LEVEL} ]; then LEVEL='0'; fi
if [ ! ${SUFFIX_ORI} ]; then SUFFIX_ORI='fastq'; fi
echo "ARGS: INPUTPATH: ${INPATH}; OUTPATH_ROOT: ${OUTPATH_ROOT}; UNIQ_ID: ${UNIQ_ID}; LENGTH: ${LENGTH}; QUAL: ${QUAL}; DEPTH: ${DEPTH}; SUFFIX_ORI: ${SUFFIX_ORI}; MULTIPRO_NUMBER: ${MULTIPRO_NUM}; DELETE_ORIGINAL: ${DELETE_ORIGINAL}; LEVEL: ${LEVEL}"
# exit 1
sjmname=`pwd`/${UNIQ_ID}.sjm
echo ${sjmname}

# basic settings
# GRCh37_75_genome="/public/home/Songlab/songyl/database/genome/human/GRCh37_release75/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
REF_FA=/public/home/Songlab/xial/Work/database/hg19/hg19_softmasked.fa
GTF_FILE=/public/home/Songlab/xial//Work/database/hg19/gencode.v27lift37.annotation.sorted.gtf.gz
ALU_BED=/public/home/Songlab/songyl/database/hg19_Alu.bed
NONALU_BED=/public/home/Songlab/songyl/database/Repeat/human/human_nonAlu.bed
EDITING_SITE=/public/home/Songlab/songyl/database/editing_sites/Merged/human/human_all_merged_editing-sites_add-strand
SIMPLE_REPEAT=/public/home/Songlab/songyl/database/hg19_SimpleRepeat.txt

if [[ ! -d ${OUTPATH_ROOT} ]]; then
    mkdir ${OUTPATH_ROOT}
fi
if [[ -f ${sjmname} ]]; then rm ${sjmname}; fi 
    touch ${sjmname}
if [[ ${DELETE_ORIGINAL} == '1' ]]; then
    if [[ -d ${OUTPATH_ROOT}/logs ]]; then rm -r ${OUTPATH_ROOT}/logs; fi
    mkdir -p ${OUTPATH_ROOT}/logs 
else
    if [[ ! -d ${OUTPATH_ROOT}/logs ]]; then mkdir ${OUTPATH_ROOT}/logs; fi
fi
if [[ ! -f ${OUTPATH_ROOT}/logs/${UNIQ_ID}.RUNNING_INFO.txt ]]; then
    touch ${OUTPATH_ROOT}/logs/${UNIQ_ID}.RUNNING_INFO.txt
fi
echo -e "\n##################################################" >> ${OUTPATH_ROOT}/logs/${UNIQ_ID}.RUNNING_INFO.txt
echo "START Running time: `/usr/bin/date`" >> ${OUTPATH_ROOT}/logs/${UNIQ_ID}.RUNNING_INFO.txt
echo -e "Paramaters:\n INPATH=${INPATH};\n OUTPATH_ROOT=${OUTPATH_ROOT};\n UNIQ_ID=${UNIQ_ID};\n QUALITY=${QUAL};\n LENGTH=${LENGTH};\n DEPTH=${DEPTH}" >> ${OUTPATH_ROOT}/logs/${UNIQ_ID}.RUNNING_INFO.txt

echo -e "job_begin\n name ${UNIQ_ID}_START\n memory 4G\n directory ${OUTPATH_ROOT}\n cmd_begin" >> ${sjmname}
echo -e "    /usr/bin/date;" >> ${sjmname}
echo -e "cmd_end\njob_end\n" >> ${sjmname}

# 00: post-run basecalling with the lastest version of Guppy after live basecalling using an older version of Guppy in MinKNOW.
# The process was based on Zhong et al., Nat Commun, 2023
if [[ ${DELETE_ORIGINAL} == '1' ]]; then
    if [[ -d ${OUTPATH_ROOT}/00_Nanoplot ]]; then rm -r ${OUTPATH_ROOT}/00_Nanoplot; fi
    mkdir -p ${OUTPATH_ROOT}/00_Nanoplot
else
    if [[ ! -d ${OUTpath} ]]; then mkdir -p ${OUTPATH_ROOT}/00_Nanoplot; fi
fi
for FASTQ in ${INPATH}/*.${SUFFIX_ORI}*;
do
    uniq_id=${UNIQ_ID}_`basename ${FASTQ} | sed "s/.${SUFFIX_ORI}//g"`
    echo -e "job_begin\n name ${uniq_id}_00Nanoplot\n memory 400G\n directory ${OUTPATH_ROOT}\n cmd_begin" >> ${sjmname}
    echo -e "   /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoPlot --fastq ${FASTQ} --loglength -o ${OUTPATH_ROOT}/00_Nanoplot --format pdf -p ${uniq_id}. --plots dot hex -t 8 > ${OUTPATH_ROOT}/00_Nanoplot/${uniq_id}_nanoplot.log 2>&1;" >> ${sjmname}
    echo -e "cmd_end\njob_end\n" >> ${sjmname}
done
# echo -e "job_begin\n name ${UNIQ_ID}_00Nanoplot\n memory 400G\n directory ${OUTPATH_ROOT}\n cmd_begin" >> ${sjmname}
# echo -e "   /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoPlot --fastq ${INPATH}/*.fastq --loglength -o ${OUTPATH_ROOT}/00_Nanoplot --format pdf -p ${UNIQ_ID}_nanoplot. --plots dot hex -t 8 > ${OUTPATH_ROOT}/00_Nanoplot/nanoplot.log 2>&1;" >> ${sjmname}
# echo -e "cmd_end\njob_end\n" >> ${sjmname}
if [[ -f ${OUTPATH_ROOT}/00_Nanoplot/${UNIQ_ID}_00Nanoplot.sh ]]; then rm ${OUTPATH_ROOT}/00_Nanoplot/${UNIQ_ID}_00Nanoplot.sh; fi
touch ${OUTPATH_ROOT}/00_Nanoplot/${UNIQ_ID}_00Nanoplot.sh
num_x=0
for FASTQ in ${INPATH}/*.${SUFFIX_ORI}*;
do
    uniq_id=${UNIQ_ID}_`basename ${FASTQ} | sed "s/.${SUFFIX_ORI}//g"`
    echo -e "/public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoPlot --fastq ${FASTQ} --loglength -o ${OUTPATH_ROOT}/00_Nanoplot --format pdf -p ${uniq_id}. --plots dot hex -t 8 > ${OUTPATH_ROOT}/00_Nanoplot/${uniq_id}_nanoplot.log 2>&1 &;" >> ${OUTPATH_ROOT}/00_Nanoplot/${UNIQ_ID}_00Nanoplot.sh
    num_x=`expr ${num_x} + 1`
    let xxx=${num_x}%${MULTIPRO_NUM}
    if [[ ${xxx} -eq 0 ]]; then 
        echo -e "wait" >> ${OUTPATH_ROOT}/00_Nanoplot/${UNIQ_ID}_00Nanoplot.sh
    fi
done
echo -e "/public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoPlot --fastq ${INPATH}/*.${SUFFIX_ORI} --loglength -o ${OUTPATH_ROOT}/00_Nanoplot --format pdf -p ${UNIQ_ID}_nanoplot. --plots dot hex -t 8 > ${OUTPATH_ROOT}/00_Nanoplot/nanoplot.log 2>&1;" >> ${OUTPATH_ROOT}/00_Nanoplot/${UNIQ_ID}_00Nanoplot.sh

#! not to run porechop, because many researches didn't conduct this
# # 01: porechop
# OUTpath=${OUTPATH_ROOT}/01_porechop
# if [[ ${DELETE_ORIGINAL} == '1' ]]; then
#     if [[ -d ${OUTpath} ]]; then rm -r ${OUTpath}; fi
#     mkdir -p ${OUTpath}
# else
#     if [[ ! -d ${OUTpath} ]]; then mkdir -p ${OUTpath}; fi
# fi
# for FASTQ in ${INPATH}/*.fastq*;
# do
#     uniq_id=${UNIQ_ID}_`basename ${FASTQ%.*} | sed 's/.fastq//g'`
#     echo -e "job_begin\n name ${uniq_id}_01porechop\n memory 400G\n directory ${OUTPATH_ROOT}\n cmd_begin" >> ${sjmname}
#     echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/porechop -i ${FASTQ} -o ${OUTpath}/${uniq_id}.trim.fastq.gz --discard_middle --format fastq.gz > ${OUTpath}/${uniq_id}.porechop.log 2>&1;" >> ${sjmname}
#     echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoPlot --fastq ${OUTpath}/${uniq_id}.trim.fastq.gz --loglength -o ${OUTpath}/nanoplot --format pdf -p ${uniq_id}. --plots dot hex -t 8 > ${OUTpath}/${uniq_id}_nanoplot.log 2>&1;" >> ${sjmname}
#     echo -e "cmd_end\njob_end\n" >> ${sjmname}
# done
# # echo -e "job_begin\n name ${UNIQ_ID}_01nanoplot\n memory 200G\n directory ${OUTpath}\n cmd_begin" >> ${sjmname}
# # echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoPlot --fastq ${OUTpath}/${UNIQ_ID}*.trim.fastq.gz --loglength -o ${OUTpath}/nanoplot --format pdf -p ${UNIQ_ID}_nanoplot. --plots dot hex -t 8 > ${OUTpath}/nanoplot.log 2>&1;" >> ${sjmname}
# # echo -e "cmd_end\njob_end\n" >> ${sjmname}
# if [[ -f ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh ]]; then rm ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh; fi
# touch ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh
# echo -e "porechop_func() {" >> ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh
# echo -e "    FASTQ=\$1" >> ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh
# echo -e "    OUTpath=\$2" >> ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh
# echo -e "    uniq_id=\$3" >> ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh
# echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/porechop -i \${FASTQ} -o \${OUTpath}/\${uniq_id}.trim.fastq.gz --discard_middle --format fastq.gz > \${OUTpath}/\${uniq_id}.porechop.log 2>&1" >> ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh
# echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoPlot --fastq \${OUTpath}/\${uniq_id}.trim.fastq.gz --loglength -o \${OUTpath}/nanoplot --format pdf -p \${uniq_id}. --plots dot hex -t 8 > \${OUTpath}/\${uniq_id}_nanoplot.log 2>&1" >> ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh
# echo -e "}" >> ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh
# num_x=0
# for FASTQ in ${INPATH}/*.fastq*;
# do
#     uniq_id=${UNIQ_ID}_`basename ${FASTQ%.*} | sed 's/.fastq//g'`
#     echo -e "porechop_func ${FASTQ} ${OUTpath} ${uniq_id} &" >> ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh
#     num_x=`expr ${num_x} + 1`
#     let xxx=${num_x}%${MULTIPRO_NUM}
#     if [[ ${xxx} -eq 0 ]]; then 
#         echo -e "wait" >> ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh
#     fi
# done
# echo -e "wait" >> ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh
# echo -e "/public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoPlot --fastq ${OUTpath}/${UNIQ_ID}*.trim.fastq.gz --loglength -o ${OUTpath}/nanoplot --format pdf -p ${UNIQ_ID}_nanoplot. --plots dot hex -t 8 > ${OUTpath}/nanoplot.log 2>&1;" >> ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh

# 02: Nanofilt
OUTpath=${OUTPATH_ROOT}/02_Nanofilt
if [[ ${DELETE_ORIGINAL} == '1' ]]; then
    if [[ -d ${OUTpath} ]]; then rm -r ${OUTpath}; fi
    mkdir -p ${OUTpath}
else
    if [[ ! -d ${OUTpath} ]]; then mkdir -p ${OUTpath}; fi
fi
for FASTQ in ${INPATH}/*.${SUFFIX_ORI}*;
do
    uniq_id=${UNIQ_ID}_`basename ${FASTQ} | sed "s/.${SUFFIX_ORI}//g"`
    echo -e "job_begin\n name ${uniq_id}_02Nanofilt\n memory 200G\n directory ${OUTPATH_ROOT}\n cmd_begin" >> ${sjmname}
    if [[ ${LENGTH} -eq 0 ]]; then
        if [[ ${SUFFIX_ORI} =~ 'gz' ]]; then
            echo -e "    gunzip -c ${FASTQ} | /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoFilt -q ${QUAL} --headcrop 50 | gzip > ${OUTpath}/${uniq_id}.trim.nanofilt.fastq.gz;" >> ${sjmname}
        else
            echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoFilt -q ${QUAL} --headcrop 50 ${FASTQ} | gzip > ${OUTpath}/${uniq_id}.trim.nanofilt.fastq.gz;" >> ${sjmname}
        fi
    else
        if [[ ${SUFFIX_ORI} =~ 'gz' ]]; then
            echo -e "    gunzip -c ${FASTQ} | /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoFilt -l ${LENGTH} -q ${QUAL} --headcrop 50 | gzip > ${OUTpath}/${uniq_id}.trim.nanofilt.fastq.gz;" >> ${sjmname}
        else
            echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoFilt -q ${QUAL} --headcrop 50 ${FASTQ} | gzip > ${OUTpath}/${uniq_id}.trim.nanofilt.fastq.gz;" >> ${sjmname}
        fi
    fi
    echo -e "    # /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/nanopolish index -s ${SEQUENCING_SUMMARY} -d ${FAST5_FILE} ${OUTpath}/02_Nanofilt/${uniq_id}.trim.nanofilt.fastq.gz > ${OUTpath}/${uniq_id}.nanopolish_index.log 2>&1" >> ${sjmname}
    echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoPlot --fastq ${OUTpath}/${uniq_id}.trim.nanofilt.fastq.gz --loglength -o ${OUTpath}/nanoplot --format pdf -p ${uniq_id}. --plots dot hex -t 8 > ${OUTpath}/${uniq_id}.nanoplot.log 2>&1;" >> ${sjmname}
    echo -e "cmd_end\njob_end\n" >> ${sjmname}
done
# echo -e "job_begin\n name ${UNIQ_ID}_02nanoplot\n memory 200G\n directory ${OUTpath}\n cmd_begin" >> ${sjmname}
# echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoPlot --fastq ${OUTpath}/${UNIQ_ID}*.trim.nanofilt.fastq.gz --loglength -o ${OUTpath}/nanoplot --format pdf -p ${UNIQ_ID}_nanoplot. --plots dot hex -t 8 > ${OUTpath}/nanoplot.log 2>&1;" >> ${sjmname}
# echo -e "cmd_end\njob_end\n" >> ${sjmname}
if [[ -f ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh ]]; then rm ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh; fi
touch ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
echo -e "nanofilt_func(){" >> ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
echo -e "    LENGTH=\$1" >> ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
echo -e "    OUTPATH_ROOT=\$2" >> ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
echo -e "    uniq_id=\$3" >> ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
echo -e "    QUAL=\$4" >> ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
echo -e "    OUTpath=\$5" >> ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
echo -e "    if [[ ${LENGTH} -eq 0 ]]; then
        if [[ ${SUFFIX_ORI} =~ 'gz' ]]; then
            gunzip -c ${FASTQ} | /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoFilt -q ${QUAL} --headcrop 50 | gzip > ${OUTpath}/${uniq_id}.trim.nanofilt.fastq.gz
        else
            /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoFilt -q ${QUAL} --headcrop 50 ${FASTQ} | gzip > ${OUTpath}/${uniq_id}.trim.nanofilt.fastq.gz
        fi
    else
        if [[ ${SUFFIX_ORI} =~ 'gz' ]]; then
            gunzip -c ${FASTQ} | /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoFilt -l ${LENGTH} -q ${QUAL} --headcrop 50 | gzip > ${OUTpath}/${uniq_id}.trim.nanofilt.fastq.gz
        else
            /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoFilt -q ${QUAL} --headcrop 50 ${FASTQ} | gzip > ${OUTpath}/${uniq_id}.trim.nanofilt.fastq.gz
        fi
    fi" >> ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoPlot --fastq ${OUTpath}/${uniq_id}.trim.nanofilt.fastq.gz --loglength -o ${OUTpath}/nanoplot --format pdf -p ${uniq_id}. --plots dot hex -t 8 > ${OUTpath}/${uniq_id}.nanoplot.log 2>&1" >> ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
echo -e "}" >> ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
num_x=0
for FASTQ in ${INPATH}/*.${SUFFIX_ORI}*;
do
    uniq_id=${UNIQ_ID}_`basename ${FASTQ} | sed "s/.${SUFFIX_ORI}//g"`
    echo -e "nanofilt_func ${LENGTH} ${OUTPATH_ROOT} ${uniq_id} ${QUAL} ${OUTpath} &" >> ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
    num_x=`expr ${num_x} + 1`
    let xxx=${num_x}%${MULTIPRO_NUM}
    if [[ ${xxx} -eq 0 ]]; then 
        echo -e "wait" >> ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
    fi
done
echo -e "wait" >> ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
echo -e "/public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoPlot --fastq ${OUTpath}/${UNIQ_ID}*.trim.nanofilt.fastq.gz --loglength -o ${OUTpath}/nanoplot --format pdf -p ${UNIQ_ID}_nanoplot. --plots dot hex -t 8 > ${OUTpath}/nanoplot.log 2>&1;" >> ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
    
# 03: minimap2 and editing sites calling
OUTpath=${OUTPATH_ROOT}/03_minimap2
if [[ ${DELETE_ORIGINAL} == '1' ]]; then
    if [[ -d ${OUTpath} ]]; then rm -r ${OUTpath}; fi
    mkdir -p ${OUTpath}
else
    if [[ ! -d ${OUTpath} ]]; then mkdir -p ${OUTpath}; fi
fi
for FASTQ in ${INPATH}/*.${SUFFIX_ORI}*;
do
    uniq_id=${UNIQ_ID}_`basename ${FASTQ} | sed "s/.${SUFFIX_ORI}//g"`
    echo -e "job_begin\n name ${uniq_id}_03minimap2\n memory 200G\n directory ${OUTPATH_ROOT}\n cmd_begin" >> ${sjmname}
    echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/minimap2 --secondary=no --cs -ax splice -uf -k14 ${REF_FA} ${OUTPATH_ROOT}/02_Nanofilt/${uniq_id}.trim.nanofilt.fastq.gz | /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/samtools view -q 30 -b -@ 4 -F 2052 | /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/samtools sort -@ 4 -O BAM -o ${OUTpath}/${uniq_id}.bam > ${OUTpath}/${uniq_id}.log 2>&1;"  >> ${sjmname} # parameter '--cs' and '-F 2052' were followed by Liu et al., Genome Biol, 2023; parameter '-ax splice -uf -k14' was followed by Workman et al., Nat Methods, 2019
    echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/samtools index ${OUTpath}/${uniq_id}.bam;" >> ${sjmname}
    echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoPlot --bam ${OUTpath}/${uniq_id}.bam --loglength -o ${OUTpath}/nanoplot --format pdf -p ${uniq_id}. --plots dot hex -t 8 > ${OUTpath}/${uniq_id}.nanoplot.log 2>&1;" >> ${sjmname}
    echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/samtools flagstat ${OUTpath}/${uniq_id}.bam > ${OUTpath}/${uniq_id}.bam.flagstat;" >> ${sjmname}
    echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/python /public/home/Songlab/xial/Work/software/L-GIREMI/bin/correct_read_strand -b ${OUTpath}/${uniq_id}.bam -o ${OUTpath}/${uniq_id}.bam -t 10 --annotation_gtf ${GTF_FILE} --genome_fasta ${REF_FA} > ${OUTpath}/${uniq_id}.correctStrand.log 2>&1;" >> ${sjmname}
    if [[ ${LEVEL} == '1' ]]; then
        echo -e "    awk '{print \$1\"\t\"\$2-1\"\t\"\$2}' ${EDITING_SITE} > ${OUTpath}/${uniq_id}_target_site.bed; " >> ${sjmname}
        echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/samtools mpileup --output-BP -A -B -d 10000 -f ${REF_FA} -q 20 -Q ${QUAL} -l ${OUTpath}/${uniq_id}_target_site.bed ${OUTpath}/${uniq_id}.bam > ${OUTpath}/${uniq_id}.mpileup;" >> ${sjmname}
        echo -e "    rm ${OUTpath}/${uniq_id}_target_site.bed;" >> ${sjmname}
        echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/python /public/home/Songlab/xial/scripts/pipeline/A2I_editing_pipe/A2I_pipe/mpileup_parser.py -i ${OUTpath}/${uniq_id}.mpileup -o ${OUTpath}/${uniq_id}.edit_level.txt --level ${EDITING_SITE} -p 10 > ${OUTpath}/${uniq_id}.mpileup_parser.log 2>&1;" >> ${sjmname}
    fi
    # For using /public/usr/bin/bcftools to conduct variants calling, Karst et al., Nat Methods, 2021 and Hall et al., Lancet Microbe, 2022 used this tools, but were both for target Nanopore sequencing or microbiome Nanopore sequencing. 
    # Medaka, another Nanopore-recommended software for variants calling, was developed only for Nanopore sequencing, some of tools require detailed version of sequencer.
    # Nanopolish, another widely-used Nanopore-specific software for data analysis including variants calling, was only developed for Nanopore, and some of tools require indexed FASTA file as input, which needs the SEQUENCING_SUMMARY file and FAST5 file that are not always provided by online datasets or PacBio sequencer.
    # echo -e "    /public/usr/bin//public/usr/bin/bcftools mpileup -B -q 30 -Q ${QUAL} -I -d 10000 -f ${REF_FA} ${OUTpath}/${uniq_id}.bam | /public/usr/bin//public/usr/bin/bcftools call -P0.01 --ploidy 1 -vm -Oz -o ${OUTpath}/${uniq_id}.vcf.gz > ${OUTpath}/${uniq_id}./public/usr/bin/bcftools.log 2>&1;" >> ${sjmname} # -Q7 may impacts the DP4 parameter. Parameter '-I', '-Q', '--ploidy' and '-m' were used according to Hall et al., Lancet Microbe, 2022. '-I' may hide some homo-edited sites?
    # echo -e "    /public/usr/bin//public/usr/bin/bcftools index ${OUTpath}/${uniq_id}.vcf.gz;" >> ${sjmname}
    echo -e "cmd_end\njob_end\n" >> ${sjmname}
done
# echo -e "job_begin\n name ${UNIQ_ID}_03nanoplot\n memory 200G\n directory ${OUTpath}\n cmd_begin" >> ${sjmname}
# echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoPlot --bam ${OUTpath}/${UNIQ_ID}*.bam --loglength -o ${OUTpath}/nanoplot --format pdf -p ${UNIQ_ID}_nanoplot. --plots dot hex -t 8 > ${OUTpath}/nanoplot.log 2>&1;" >> ${sjmname}
# echo -e "cmd_end\njob_end\n" >> ${sjmname}
if [[ -f ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03_minimap2.sh ]]; then rm ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03_minimap2.sh; fi
touch ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03_minimap2.sh
echo -e "minimap_func() {" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03_minimap2.sh
echo -e "    uniq_id=\$1" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03_minimap2.sh
echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/minimap2 --secondary=no --cs -ax splice -uf -k14 ${REF_FA} ${OUTPATH_ROOT}/02_Nanofilt/\${uniq_id}.trim.nanofilt.fastq.gz | /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/samtools view -q 30 -b -@ 4 -F 2052 | /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/samtools sort -@ 4 -O BAM -o ${OUTpath}/\${uniq_id}.bam > ${OUTpath}/\${uniq_id}.log 2>&1;"  >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03_minimap2.sh # parameter '--cs' and '-F 2052' were followed by Liu et al., Genome Biol, 2023; parameter '-ax splice -uf -k14' was followed by Workman et al., Nat Methods, 2019
echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/samtools index ${OUTpath}/\${uniq_id}.bam;" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03_minimap2.sh
echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoPlot --bam ${OUTpath}/\${uniq_id}.bam --loglength -o ${OUTpath}/nanoplot --format pdf -p \${uniq_id}. --plots dot hex -t 8 > ${OUTpath}/\${uniq_id}.nanoplot.log 2>&1;" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03_minimap2.sh
echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/samtools flagstat ${OUTpath}/\${uniq_id}.bam > ${OUTpath}/\${uniq_id}.bam.flagstat;" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03_minimap2.sh
echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/python /public/home/Songlab/xial/Work/software/L-GIREMI/bin/correct_read_strand -b ${OUTpath}/\${uniq_id}.bam -o ${OUTpath}/\${uniq_id}.bam -t 10 --annotation_gtf ${GTF_FILE} --genome_fasta ${REF_FA} > ${OUTpath}/\${uniq_id}.correctStrand.log 2>&1;" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03_minimap2.sh
if [[ ${LEVEL} == '1' ]]; then
    echo -e "    awk '{print \$1\"\t\"\$2-1\"\t\"\$2}' ${EDITING_SITE} > ${OUTpath}/\${uniq_id}_target_site.bed" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03_minimap2.sh
    echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/samtools mpileup --output-BP -A -B -d 10000 -f ${REF_FA} -q 20 -Q ${QUAL} -l ${OUTpath}/\${uniq_id}_target_site.bed ${OUTpath}/\${uniq_id}.bam > ${OUTpath}/\${uniq_id}.mpileup" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03_minimap2.sh
    echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/python /public/home/Songlab/xial/scripts/pipeline/A2I_editing_pipe/A2I_pipe/mpileup_parser.py -i ${OUTpath}/\${uniq_id}.mpileup -o ${OUTpath}/\${uniq_id}.edit_level.txt --level ${EDITING_SITE} -p 5 > ${OUTpath}/\${uniq_id}.mpileup_parser.log 2>&1" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03_minimap2.sh
    echo -e "    rm ${OUTpath}/\${uniq_id}_target_site.bed;" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03_minimap2.sh
fi
echo -e "}" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03_minimap2.sh
for FASTQ in ${INPATH}/*.${SUFFIX_ORI}*;
do
    uniq_id=${UNIQ_ID}_`basename ${FASTQ} | sed "s/.${SUFFIX_ORI}//g"`
    echo -e "minimap_func ${uniq_id} &" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03_minimap2.sh
    # echo -e "/public/usr/bin//public/usr/bin/bcftools mpileup -B -q ${QUAL} -Q 0 -I -d 10000 -f ${REF_FA} ${OUTpath}/${uniq_id}.bam | /public/usr/bin//public/usr/bin/bcftools call -P0.01 --ploidy 1 -vm -Oz -o ${OUTpath}/${uniq_id}.vcf.gz > ${OUTpath}/${uniq_id}./public/usr/bin/bcftools.log 2>&1;" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03_minimap2.sh # -Q7 may impacts the DP4 parameter. Parameter '-I', '-Q', '--ploidy' and '-m' were used according to Hall et al., Lancet Microbe, 2022. '-I' may hide some homo-edited sites?
    # nohup samtools mpileup -A -B -d 10000 -f /public/home/Songlab/songyl/database/ref_genomes/Hg19/hg19_softmasked.fa -q 30 -Q 7 -o SRP275847_nanopore_SRR12389273.mpileup ../SRP275847_nanopore_SRR12389273.bam > SRP275847_nanopore_SRR12389273.mpileup.log 2>&1 & ## /public/usr/bin/bcftools mpileupu + call skipped many sites edited in cluster, while samtools mpileup stored such sites
    # echo -e "/public/usr/bin//public/usr/bin/bcftools index ${OUTpath}/${uniq_id}.vcf.gz;" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03_minimap2.sh
done
echo -e "wait" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03_minimap2.sh
# if [[ ${LEVEL} == '1' ]]; then
#     echo -e "/usr/bin/perl /public/home/Songlab/xial/scripts/pipeline_songyl/call_levels/editing_site_call_from_RNAseq_ADAR_20170512.pl --in ${OUTPATH_ROOT}/03_minimap2 --out ${OUTPATH_ROOT}/03_minimap2/callKnownSitesLevels > mpileup_${UNIQ_ID}.log 2>&1;" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03_minimap2.sh
# fi
echo -e "/public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoPlot --bam ${OUTpath}/${UNIQ_ID}*.bam --loglength -o ${OUTpath}/nanoplot --format pdf -p ${UNIQ_ID}_nanoplot. --plots dot hex -t 8 > ${OUTpath}/nanoplot.log 2>&1;" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03_minimap2.sh
    

# 04 call sites
if [[ ${LEVEL} == '0' ]]; then
    # 04: call sites separated according to chromosome
    OUTpath=${OUTPATH_ROOT}/04_callSite
    snp_db="/public/home/Songlab/xial/Work/database/hg19/SNP_all_dbSNP_1000Genomes_UWash.txt"
    if [[ ${DELETE_ORIGINAL} == '1' ]]; then
        if [[ -d ${OUTpath} ]]; then rm -r ${OUTpath}; fi
        mkdir -p ${OUTpath}
    else
        if [[ ! -d ${OUTpath} ]]; then mkdir -p ${OUTpath}; fi
    fi
    chr_list=(`seq 1 22` 'X' 'Y')
    for FASTQ in ${INPATH}/*.${SUFFIX_ORI}*;
    do
        uniq_id=${UNIQ_ID}_`basename ${FASTQ} | sed "s/.${SUFFIX_ORI}//g"`
        echo -e "job_begin\n name ${uniq_id}_04_bcftools_call\n memory 200G\n directory ${OUTPATH_ROOT}\n cmd_begin" >> ${sjmname}
        echo -e "    /public/usr/bin/bcftools mpileup --threads 8 -Ou -f ~/Work/database/hg19/hg19_softmasked.fa ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_`basename ${FASTQ%.*} | sed "s/.${SUFFIX_ORI}//g"`.bam | /public/usr/bin/bcftools call --ploidy 1 --threads 8 -mv -Oz -o ${OUTpath}/${uniq_id}.vcf.gz" >> ${sjmnames}
        echo -e "    /public/usr/bin/bedtools intersect -a ${OUTpath}/${uniq_id}.vcf.gz -b /public/home/Songlab/songyl/database/hg19_Alu.bed -wa -header > ${OUTpath}/${uniq_id}.ALU.vcf" >> ${sjmname}
        echo -e "    /public/usr/bin/bedtools intersect -v -a ${OUTpath}/${uniq_id}.ALU.vcf -b ../../04_callSite/model_construction_test/snp_db.bed -wa -header >  ${OUTpath}/${uniq_id}.ALU.remove_snp.vcf" >> ${sjmname}
        echo -e "    /public/usr/bin/bcftools filter --threads 8 -e 'MQ < 10.0 || QUAL < 10' -Oz -o ${OUTpath}/${uniq_id}.ALU.remove_snp.filt.vcf.gz ${OUTpath}/${uniq_id}.ALU.remove_snp.vcf" >> ${sjmname}
        echo -e "    /public/usr/bin/bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%DP\t%VDB\t%RPBZ\t%MQBZ\t%NMBZ\t%SCBZ\t%MQ0F\t%DP4\t%MQ\n' -o ${OUTpath}/${uniq_id}.ALU.cand_site.remove_snp.filt.vcfTable.tmp.txt ${OUTpath}/${uniq_id}.ALU.remove_snp.filt.vcf.gz" >> ${sjmname}
        echo -e "    awk -v OFS='\t' '{if(length($3)==1 && length($4)==1){split($12,str_x,","); print $0,str_x[3]+str_x[4]}}' ${OUTpath}/${uniq_id}.ALU.cand_site.remove_snp.filt.vcfTable.tmp.txt | sed '1i CHROM\tPOS\tREF\tALT\tDP\tVDB\tRPBZ\tMQBZ\tNMBZ\tSCBZ\tMQ0F\tDP4\tMQ\tMN' > ${OUTpath}/${uniq_id}.ALU.cand_site.remove_snp.filt.vcfTable.txt" >> ${sjmname}
        echo -e "cmd_end\njob_end\n" >> ${sjmname}
        for key in ${chr_list[@]}; do
            uniq_id=${UNIQ_ID}_`basename ${FASTQ} | sed "s/.${SUFFIX_ORI}//g"`_chr${key}
            echo -e "job_begin\n name ${uniq_id}_04_mpileup_call\n memory 200G\n directory ${OUTPATH_ROOT}\n cmd_begin" >> ${sjmname}
            echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/samtools view -q 30 -hb -@ 4 -F 2052 ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_`basename ${FASTQ%.*} | sed "s/.${SUFFIX_ORI}//g"`.bam chr$key | /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/samtools sort -@ 4 -o ${OUTpath}/${uniq_id}.bam;" >> ${sjmname}
            echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/samtools index ${OUTpath}/${uniq_id}.bam;" >> ${sjmname}
            echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/samtools mpileup --output-BP -A -B -d 10000 -f ${REF_FA} -q 30 -Q ${QUAL} -l ${ALU_BED} -o ${OUTpath}/${uniq_id}.ALU.mpileup ${OUTpath}/${uniq_id}.bam;" >> ${sjmname}
            echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/python /public/home/Songlab/xial/scripts/pipeline/A2I_editing_pipe/A2I_pipe/mpileup_parser.py -i ${OUTpath}/${uniq_id}.ALU.mpileup -o ${OUTpath}/${uniq_id}.ALU.candsite.txt -R ${REF_FA} -Q ${QUAL} -sr ${SIMPLE_REPEAT} > ${OUTpath}/${uniq_id}.ALU.mpileup_parser.log 2>&1;" >> ${sjmname}
            echo -e "cmd_end\njob_end\n" >> ${sjmname}
        done
    done
    if [[ -f ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh ]]; then rm ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh; fi
    touch ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
    for FASTQ in ${INPATH}/*.${SUFFIX_ORI}*;
    do
        uniq_id=${UNIQ_ID}_`basename ${FASTQ} | sed "s/.${SUFFIX_ORI}//g"`
        echo -e "/public/usr/bin/bcftools mpileup --threads 8 -Ou -f ~/Work/database/hg19/hg19_softmasked.fa ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_\`basename ${FASTQ%.*} | sed "s/.${SUFFIX_ORI}//g"\`.bam | /public/usr/bin/bcftools call --ploidy 1 --threads 8 -mv -Oz -o ${OUTpath}/${uniq_id}.vcf.gz" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
        echo -e "/public/usr/bin/bedtools intersect -a ${OUTpath}/${uniq_id}.vcf.gz -b /public/home/Songlab/songyl/database/hg19_Alu.bed -wa -header > ${OUTpath}/${uniq_id}.ALU.vcf" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
        echo -e "/public/usr/bin/bedtools intersect -v -a ${OUTpath}/${uniq_id}.ALU.vcf -b ../../04_callSite/model_construction_test/snp_db.bed -wa -header >  ${OUTpath}/${uniq_id}.ALU.remove_snp.vcf" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
        echo -e "/public/usr/bin/bcftools filter --threads 8 -e 'MQ < 10.0 || QUAL < 10' -Oz -o ${OUTpath}/${uniq_id}.ALU.remove_snp.filt.vcf.gz ${OUTpath}/${uniq_id}.ALU.remove_snp.vcf" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
        echo -e "/public/usr/bin/bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%DP\t%VDB\t%RPBZ\t%MQBZ\t%NMBZ\t%SCBZ\t%MQ0F\t%DP4\t%MQ\n' -o ${OUTpath}/${uniq_id}.ALU.cand_site.remove_snp.filt.vcfTable.tmp.txt ${OUTpath}/${uniq_id}.ALU.remove_snp.filt.vcf.gz" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
        echo -e "awk -v OFS='\t' '{if(length($3)==1 && length($4)==1){split($12,str_x,","); print $0,str_x[3]+str_x[4]}}' ${OUTpath}/${uniq_id}.ALU.cand_site.remove_snp.filt.vcfTable.tmp.txt | sed '1i CHROM\tPOS\tREF\tALT\tDP\tVDB\tRPBZ\tMQBZ\tNMBZ\tSCBZ\tMQ0F\tDP4\tMQ\tMN' > ${OUTpath}/${uniq_id}.ALU.cand_site.remove_snp.filt.vcfTable.txt" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
        echo -e "\n" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
    done
    echo -e "\n\n" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
    echo -e "chr_call_site() {" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
    echo -e "    FASTQ=\$1" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
    echo -e "    uniq_id=\$2" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
    echo -e "    key=\$3" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
    echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/samtools view -q 30 -hb -@ 4 -F 2052 ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_\`basename \${FASTQ%.*} | sed 's/.fastq//g'\`.bam chr\${key} | /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/samtools sort -@ 4 -o ${OUTpath}/\${uniq_id}.bam;" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
    echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/samtools index ${OUTpath}/\${uniq_id}.bam;" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
    echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/samtools mpileup --output-BP -A -B -d 10000 -f ${REF_FA} -q 30 -Q ${QUAL} -l ${ALU_BED} -o  ${OUTpath}/\${uniq_id}.ALU.mpileup ${OUTpath}/\${uniq_id}.bam;" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
    echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/python /public/home/Songlab/xial/scripts/pipeline/A2I_editing_pipe/A2I_pipe/mpileup_parser.py -i ${OUTpath}/\${uniq_id}.ALU.mpileup -o ${OUTpath}/\${uniq_id}.ALU.candsite.txt -Q ${QUAL} -R ${REF_FA} -sr ${SIMPLE_REPEAT} > ${OUTpath}/\${uniq_id}.ALU.mpileup_parser.log 2>&1;" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
    echo -e "}" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
    for FASTQ in ${INPATH}/*.${SUFFIX_ORI}*;
    do
        for key in ${chr_list[@]}; do
            uniq_id=${UNIQ_ID}_`basename ${FASTQ} | sed "s/.${SUFFIX_ORI}//g"`_chr${key}
            echo -e "chr_call_site ${FASTQ} ${uniq_id} ${key} &" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
        done
        echo -e 'wait' >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
    done
    echo -e "wait" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh

    # cat output
    for FASTQ in ${INPATH}/*.${SUFFIX_ORI}*;
    do
        uniq_id=${UNIQ_ID}_`basename ${FASTQ} | sed "s/.${SUFFIX_ORI}//g"`
        echo -e "job_begin\n name ${uniq_id}_04_callSite_CAT\n memory 100G\n directory ${OUTPATH_ROOT}\n cmd_begin" >> ${sjmname}
        echo -e "    cat ${OUTPATH_ROOT}/04_callSite/${uniq_id}_chr*.ALU.candsite.txt | grep -v 'CHROM' | sed '1i CHROM\tPOS\tREADS_NUM\tREF\tALT\tALT_QUALITY\tEDITLEVEL\tCOVERAGE\tCOVERAGE_noSKIP\tDELETION_NUM\tBP_RATIO\tSKIP\tPUZZLE_NUM\tSIMPLE_REPEAT\tHOMOPOLYMER\tSPLICE_NT\tSEQ_ERROR_P' > ${OUTPATH_ROOT}/04_callSite/${uniq_id}.ALU.candsite.txt;" >> ${sjmname}
        echo -e "    awk 'ARGIND==1{a[\$1\"_\"\$2]=\$1\"_\"\$2} ARGIND==2{if(!(\$1\"_\"\$2 in a)){print \$0}}' ${snp_db} ${OUTPATH_ROOT}/04_callSite/${uniq_id}.ALU.candsite.txt > ${OUTPATH_ROOT}/04_callSite/${uniq_id}.ALU.de_snpdb.candsite.txt;" >> ${sjmname}
        echo -e "    site_num=\`sed '1d' ${OUTPATH_ROOT}/04_callSite/${uniq_id}.ALU.de_snpdb.candsite.txt | wc -l\`;" >> ${sjmname}
        echo -e "    sed '1d' ${OUTPATH_ROOT}/04_callSite/${uniq_id}.ALU.de_snpdb.candsite.txt | awk '{print toupper(\$4)\">\"toupper(\$5)}' | sort | uniq -c | awk -v OFS='\t' -v site_num=\${site_num} '{print \$2,\$1,\$1/site_num}' > ${OUTPATH_ROOT}/04_callSite/${uniq_id}.ALU.de_snpdb.edit_type.txt;" >> ${sjmname}
        echo -e " cmd_end\njob_end\n" >> ${sjmname}
    done
    echo -e "CAT_SITE() {" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
    echo -e "    uniq_id=\$1" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
    echo -e "    cat ${OUTPATH_ROOT}/04_callSite/\${uniq_id}_chr*.ALU.candsite.txt | grep -v 'CHROM' | sed '1i CHROM\tPOS\tREADS_NUM\tREF\tALT\tALT_QUALITY\tEDITLEVEL\tCOVERAGE\tCOVERAGE_noSKIP\tDELETION_NUM\tBP_RATIO\tSKIP\tPUZZLE_NUM\tSIMPLE_REPEAT\tHOMOPOLYMER\tSPLICE_NT\tSEQ_ERROR_P' > ${OUTPATH_ROOT}/04_callSite/\${uniq_id}.ALU.candsite.txt;" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
    echo -e "    awk 'ARGIND==1{a[\$1\"_\"\$2]=\$1\"_\"\$2} ARGIND==2{if(!(\$1\"_\"\$2 in a)){print \$0}}' ${snp_db} ${OUTPATH_ROOT}/04_callSite/\${uniq_id}.ALU.candsite.txt > ${OUTPATH_ROOT}/04_callSite/\${uniq_id}.ALU.de_snpdb.candsite.txt;" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
    echo -e "    site_num=\`sed '1d' ${OUTPATH_ROOT}/04_callSite/\${uniq_id}.ALU.de_snpdb.candsite.txt | wc -l\`;" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
    echo -e "    sed '1d' ${OUTPATH_ROOT}/04_callSite/\${uniq_id}.ALU.de_snpdb.candsite.txt | awk '{print toupper(\$4)\">\"toupper(\$5)}' | sort | uniq -c | awk -v OFS='\t' -v site_num=\${site_num} '{print \$2,\$1,\$1/site_num}' > ${OUTPATH_ROOT}/04_callSite/\${uniq_id}.ALU.de_snpdb.edit_type.txt;" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh    
    echo -e "}" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
    for FASTQ in ${INPATH}/*.${SUFFIX_ORI}*;
    do
        uniq_id=${UNIQ_ID}_`basename ${FASTQ} | sed "s/.${SUFFIX_ORI}//g"`
        echo -e "CAT_SITE ${uniq_id} &" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
    done
    echo -e "wait" >> ${OUTPATH_ROOT}/04_callSite/${UNIQ_ID}_04_callSite.sh
fi

# # 04: dividing sites into Alu, nonAlu-repeat, other
# snp_db="/public/home/Songlab/xial/Work/database/hg19/SNP_all_dbSNP_1000Genomes_UWash.txt"
# OUTpath=${OUTPATH_ROOT}/04_region
# if [[ ${DELETE_ORIGINAL} == '1' ]]; then
#     if [[ -d ${OUTpath} ]]; then rm -r ${OUTpath}; fi
#     if [[ -f ${OUTpath}/${uniq_id}.count ]]; then rm -r ${OUTpath}/${uniq_id}.count; fi
#     mkdir -p ${OUTpath}
# fi
# count_dbSNP_depth() {
#     IN_VCF=$1
#     sjmname=$2
#     OUTpath=$3
#     uniq_id=$4
#     region=$5
#     echo -e "    original=\`less ${IN_VCF} | grep -v '##' | wc -l\`;" >> ${sjmname}
#     echo -e "    /public/usr/bin//public/usr/bin/bcftools filter --threads 4 --mask-file ~/Work/database/hg19/SNP_all_dbSNP_1000Genomes_UWash.txt -s SNP ${IN_VCF} | grep -v 'SNP' | gzip -c > ${OUTpath}/${uniq_id}.${region}.dbSNP_filt.vcf.gz;" >> ${sjmname}
#     echo -e "    dbSNP_filt=\`less ${OUTpath}/${uniq_id}.${region}.dbSNP_filt.vcf.gz | grep -v '##' | wc -l\`;" >> ${sjmname}
#     echo -e "    /public/usr/bin//public/usr/bin/bcftools filter --threads 4 --exclude 'DP4[0]+DP4[1]+DP4[2]+DP4[3]<10' -s DEPTH_FAILED ${OUTpath}/${uniq_id}.${region}.dbSNP_filt.vcf.gz | grep -v 'DEPTH_FAILED' | gzip -c > ${OUTpath}/${uniq_id}.${region}.dbSNP_depth_filt.vcf.gz;" >> ${sjmname}
#     echo -e "    depth_filt=\`less ${OUTpath}/${uniq_id}.${region}.dbSNP_depth_filt.vcf.gz | grep -v '##' | wc -l\`;" >> ${sjmname}
#     echo -e "    echo -e \"dbSNP number for `basename ${IN_VCF}`: \`expr \${original} - \${dbSNP_filt}\`\" >> ${OUTpath}/${uniq_id}.count;" >> ${sjmname}
#     echo -e "    echo -e \"depth filtered number for `basename ${IN_VCF}`: \`expr \${dbSNP_filt} - \${depth_filt}\`\" >> ${OUTpath}/${uniq_id}.count;" >> ${sjmname}
# }
# for FASTQ in ${INPATH}/*.fastq;
# do
#     uniq_id=${UNIQ_ID}_`basename ${FASTQ%.*}`
#     prefix=${OUTpath}/${uniq_id}
#     echo -e "job_begin\n name ${uniq_id}_04region\n memory 200G\n directory ${OUTPATH_ROOT}\n cmd_begin" >> ${sjmname}
#     echo -e "    /public/home/Songlab/xial/miniconda3/bin/vcftools --exclude-positions ${snp_db} --gzvcf ${OUTPATH_ROOT}/03_minimap2/${uniq_id}.vcf.gz --recode --recode-INFO-all --stdout | /public/usr/bin//public/usr/bin/bcftools view -Oz -o ${prefix}.snpdb.vcf.gz > ${prefix}.ALU.snpdb.log 2>&1;" >> ${sjmname}
#     # ALU region
#     echo -e "    /public/home/Songlab/xial/miniconda3/bin/vcftools --bed ${ALU_BED} --gzvcf ${prefix}.snpdb.vcf.gz --recode-INFO-all --recode --stdout | /public/usr/bin//public/usr/bin/bcftools view -Oz -o ${prefix}.ALU.vcf.gz > ${prefix}.ALU.gz.log 2>&1;" >> ${sjmname}
#     echo -e "    /public/usr/bin//public/usr/bin/bcftools filter --threads 4 --include 'DP4[0]+DP4[1]+DP4[2]+DP4[3]>${DEPTH}' -o ${prefix}.ALU.filtered.bcf -Ob ${prefix}.ALU.vcf.gz > ${prefix}.ALU.filtered.log 2>&1;" >> ${sjmname}
#     echo -e "    /public/usr/bin//public/usr/bin/bcftools index ${prefix}.ALU.filtered.bcf;" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh 
#     echo -e "    /public/home/Songlab/xial/miniconda3/bin/vcftools --exclude-bed ${ALU_BED} --gzvcf ${prefix}.snpdb.vcf.gz --recode-INFO-all --recode --out ${prefix}.nonALU > ${prefix}.nonALU.log 2>&1;" >> ${sjmname}
#     # nonALU repeat region
#     echo -e "    /public/home/Songlab/xial/miniconda3/bin/vcftools --bed ${NONALU_BED} --vcf ${prefix}.nonALU.recode.vcf --recode-INFO-all --recode --stdout | /public/usr/bin//public/usr/bin/bcftools view -Oz -o ${prefix}.nonALU_rep.vcf.gz > ${prefix}.nonALU_rep.gz.log 2>&1;" >> ${sjmname}
#     echo -e "    /public/usr/bin//public/usr/bin/bcftools filter --threads 4 --include 'DP4[0]+DP4[1]+DP4[2]+DP4[3]>${DEPTH}' -o ${prefix}.nonALU_rep.filtered.bcf -Ob ${prefix}.nonALU_rep.vcf.gz > ${prefix}.nonALU_rep.filtered.log 2>&1;" >> ${sjmname}
#     echo -e "    /public/usr/bin//public/usr/bin/bcftools index ${prefix}.nonALU_rep.filtered.bcf;" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
#     # other
#     echo -e "    /public/home/Songlab/xial/miniconda3/bin/vcftools --exclude-bed ${NONALU_BED} --vcf ${prefix}.nonALU.recode.vcf --recode-INFO-all --recode --stdout | /public/usr/bin//public/usr/bin/bcftools view -Oz -o ${prefix}.other.vcf.gz > ${prefix}.other.gz.log 2>&1;" >> ${sjmname}
#     echo -e "    /public/usr/bin//public/usr/bin/bcftools filter --threads 4 --include 'DP4[0]+DP4[1]+DP4[2]+DP4[3]>${DEPTH}' -o ${prefix}.other.filtered.bcf -Ob ${prefix}.other.vcf.gz > ${prefix}.other.filtered.log 2>&1;" >> ${sjmname}
#     echo -e "    /public/usr/bin//public/usr/bin/bcftools index ${prefix}.other.filtered.bcf;" >> ${sjmname}
#     # calculate
#     echo -e "    ALU_num=\`/public/usr/bin//public/usr/bin/bcftools view ${prefix}.ALU.filtered.bcf | grep -v '##' | wc -l\`;" >> ${sjmname}
#     echo -e "    nonALU_rep_num=\`/public/usr/bin//public/usr/bin/bcftools view ${prefix}.nonALU_rep.filtered.bcf | grep -v '##' | wc -l\`;" >> ${sjmname}
#     echo -e "    other_num=\`/public/usr/bin//public/usr/bin/bcftools view ${prefix}.other.filtered.bcf | grep -v '##' | wc -l\`;" >> ${sjmname}
#     echo -e "    echo -e \"ALU region number: ${ALU_num}\" >> ${prefix}.count;" >> ${sjmname}
#     echo -e "    echo -e \"nonALU_rep region number: ${nonALU_rep_num}\" >> ${prefix}.count;" >> ${sjmname}
#     echo -e "    echo -e \"Other region number: ${other_num}\" >> ${prefix}.count;" >> ${sjmname}
#     echo -e "cmd_end\njob_end\n" >> ${sjmname}
# done
# if [[ -f ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh ]]; then rm ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh; fi
# touch ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
# for FASTQ in ${INPATH}/*.fastq;
# do
#     uniq_id=${UNIQ_ID}_`basename ${FASTQ%.*}`
#     echo -e "    /public/home/Songlab/xial/miniconda3/bin/vcftools --exclude-positions ${snp_db} --gzvcf ${OUTPATH_ROOT}/03_minimap2/${uniq_id}.vcf.gz --recode --recode-INFO-all --stdout | /public/usr/bin//public/usr/bin/bcftools view -Oz -o ${prefix}.snpdb.vcf.gz > ${prefix}.ALU.snpdb.log 2>&1" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
#     # ALU region
#     echo -e "    /public/home/Songlab/xial/miniconda3/bin/vcftools --bed ${ALU_BED} --gzvcf ${prefix}.snpdb.vcf.gz --recode-INFO-all --recode --stdout | /public/usr/bin//public/usr/bin/bcftools view -Oz -o ${prefix}.ALU.vcf.gz > ${prefix}.ALU.gz.log 2>&1" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
#     echo -e "    /public/usr/bin//public/usr/bin/bcftools filter --threads 4 --include 'DP4[0]+DP4[1]+DP4[2]+DP4[3]>${DEPTH}' -o ${prefix}.ALU.filtered.bcf -Ob ${prefix}.ALU.vcf.gz > ${prefix}.ALU.filtered.log 2>&1" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
#     echo -e "    /public/usr/bin//public/usr/bin/bcftools index ${prefix}.ALU.filtered.bcf" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh 
#     echo -e "    /public/home/Songlab/xial/miniconda3/bin/vcftools --exclude-bed ${ALU_BED} --gzvcf ${prefix}.snpdb.vcf.gz --recode-INFO-all --recode --out ${prefix}.nonALU > ${prefix}.nonALU.log 2>&1" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
#     # nonALU repeat region
#     echo -e "    /public/home/Songlab/xial/miniconda3/bin/vcftools --bed ${NONALU_BED} --vcf ${prefix}.nonALU.recode.vcf --recode-INFO-all --recode --stdout | /public/usr/bin//public/usr/bin/bcftools view -Oz -o ${prefix}.nonALU_rep.vcf.gz > ${prefix}.nonALU_rep.gz.log 2>&1" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
#     echo -e "    /public/usr/bin//public/usr/bin/bcftools filter --threads 4 --include 'DP4[0]+DP4[1]+DP4[2]+DP4[3]>${DEPTH}' -o ${prefix}.nonALU_rep.filtered.bcf -Ob ${prefix}.nonALU_rep.vcf.gz > ${prefix}.nonALU_rep.filtered.log 2>&1" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
#     echo -e "    /public/usr/bin//public/usr/bin/bcftools index ${prefix}.nonALU_rep.filtered.bcf" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
#     # other
#     echo -e "    /public/home/Songlab/xial/miniconda3/bin/vcftools --exclude-bed ${NONALU_BED} --vcf ${prefix}.nonALU.recode.vcf --recode-INFO-all --recode --stdout | /public/usr/bin//public/usr/bin/bcftools view -Oz -o ${prefix}.other.vcf.gz > ${prefix}.other.gz.log 2>&1" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
#     echo -e "    /public/usr/bin//public/usr/bin/bcftools filter --threads 4 --include 'DP4[0]+DP4[1]+DP4[2]+DP4[3]>${DEPTH}' -o ${prefix}.other.filtered.bcf -Ob ${prefix}.other.vcf.gz > ${prefix}.other.filtered.log 2>&1" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
#     echo -e "    /public/usr/bin//public/usr/bin/bcftools index ${prefix}.other.filtered.bcf" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
#     # calculate
#     echo -e "    ALU_num=\`/public/usr/bin//public/usr/bin/bcftools view ${prefix}.ALU.filtered.bcf | grep -v '##' | wc -l\`" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
#     echo -e "    nonALU_rep_num=\`/public/usr/bin//public/usr/bin/bcftools view ${prefix}.nonALU_rep.filtered.bcf | grep -v '##' | wc -l\`" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
#     echo -e "    other_num=\`/public/usr/bin//public/usr/bin/bcftools view ${prefix}.other.filtered.bcf | grep -v '##' | wc -l\`" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
#     echo -e "    echo -e \"ALU region number: ${ALU_num}\" >> ${prefix}.count" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
#     echo -e "    echo -e \"nonALU_rep region number: ${nonALU_rep_num}\" >> ${prefix}.count" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
#     echo -e "    echo -e \"Other region number: ${other_num}\" >> ${prefix}.count" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
# done

# # 05 candidate site filtration
# #snp_db="/public/home/Songlab/xial/Work/database/hg19/dbSNP151_SNV.site.txt,/public/home/Songlab/songyl/database/hg19_1000genomes_newSNPCalls.txt,/public/home/Songlab/songyl/database/UWash_variants.txt"
# ref_fa="/public/home/Songlab/songyl/database/ref_genomes/Hg19/hg19_softmasked.fa"
# simple_repeat="/public/home/Songlab/songyl/database/hg19_SimpleRepeat.txt"
# OUTpath=${OUTPATH_ROOT}/05_cand_sites
# if [[ ${DELETE_ORIGINAL} == '1' ]]; then
#     if [[ -d ${OUTpath} ]]; then rm -r ${OUTpath}; fi
#     mkdir -p ${OUTpath}
# fi
# for FASTQ in ${INPATH}/*.fastq;
# do
#     uniq_id=${UNIQ_ID}_`basename ${FASTQ%.*}`
#     #prefix_dir=${OUTpath}/${uniq_id}
#     echo -e "job_begin\n name ${uniq_id}_05candSites\n memory 200G\n directory ${OUTPATH_ROOT}\n cmd_begin" >> ${sjmname}
#     echo -e "    /public/home/Songlab/xial/miniconda3/bin/python /public/home/Songlab/xial/scripts/pipeline/A2I_editing_pipe/A2I_pipe/vcf2RNAvcf.py -v ${OUTPATH_ROOT}/04_region/${uniq_id}.ALU.filtered.bcf -d ${DEPTH} -r ALU -o ${OUTpath}/${uniq_id}_Alu.all_candsites > ${OUTpath}/${uniq_id}.ALU_sites.log 2>&1;" >> ${sjmname}
#     echo -e "    /public/home/Songlab/xial/miniconda3/bin/python /public/home/Songlab/xial/scripts/pipeline/A2I_editing_pipe/A2I_pipe/vcf2RNAvcf.py -v ${OUTPATH_ROOT}/04_region/${uniq_id}.nonALU_rep.filtered.bcf -d ${DEPTH} -r NONALU -o ${OUTpath}/${uniq_id}_NoAlu_rep.all_candsites -R ${REF_FA} -sr ${simple_repeat} > ${OUTpath}/${uniq_id}.nonALU_sites.log 2>&1;" >> ${sjmname}
#     echo -e "    /public/home/Songlab/xial/miniconda3/bin/python /public/home/Songlab/xial/scripts/pipeline/A2I_editing_pipe/A2I_pipe/vcf2RNAvcf.py -v ${OUTPATH_ROOT}/04_region/${uniq_id}.other.filtered.bcf -d ${DEPTH} -r NONALU -o ${OUTpath}/${uniq_id}_other.all_candsites -R ${REF_FA} -sr ${simple_repeat} > ${OUTpath}/${uniq_id}.other_sites.log 2>&1;" >> ${sjmname}
#     echo -e "cmd_end\njob_end\n" >> ${sjmname}
# done
# if [[ -f ${OUTPATH_ROOT}/05_cand_sites/${UNIQ_ID}_05candSites.sh ]]; then rm ${OUTPATH_ROOT}/05_cand_sites/${UNIQ_ID}_05candSites.sh; fi
# touch ${OUTPATH_ROOT}/05_cand_sites/${UNIQ_ID}_05candSites.sh
# for FASTQ in ${INPATH}/*.fastq;
# do
#     uniq_id=${UNIQ_ID}_`basename ${FASTQ%.*}`
#     OUTpath=${OUTPATH_ROOT}/05_cand_sites
#     #prefix_dir=${OUTpath}/${uniq_id}
#     echo -e "/public/home/Songlab/xial/miniconda3/bin/python /public/home/Songlab/xial/scripts/pipeline/A2I_editing_pipe/A2I_pipe/vcf2RNAvcf.py -v ${OUTPATH_ROOT}/04_region/${uniq_id}.ALU.filtered.bcf -d ${DEPTH} -r ALU -o ${OUTpath}/${uniq_id}_Alu.all_candsites > ${OUTpath}/${uniq_id}.ALU_sites.log 2>&1;" >> ${OUTPATH_ROOT}/05_cand_sites/${UNIQ_ID}_05candSites.sh
#     echo -e "/public/home/Songlab/xial/miniconda3/bin/python /public/home/Songlab/xial/scripts/pipeline/A2I_editing_pipe/A2I_pipe/vcf2RNAvcf.py -v ${OUTPATH_ROOT}/04_region/${uniq_id}.nonALU_rep.filtered.bcf -d ${DEPTH} -r NONALU -o ${OUTpath}/${uniq_id}_NoAlu_rep.all_candsites -R ${ref_fa} -sr ${simple_repeat} > ${OUTpath}/${uniq_id}.nonALU_sites.log 2>&1;" >> ${OUTPATH_ROOT}/05_cand_sites/${UNIQ_ID}_05candSites.sh
#     echo -e "/public/home/Songlab/xial/miniconda3/bin/python /public/home/Songlab/xial/scripts/pipeline/A2I_editing_pipe/A2I_pipe/vcf2RNAvcf.py -v ${OUTPATH_ROOT}/04_region/${uniq_id}.other.filtered.bcf -d ${DEPTH} -r NONALU -o ${OUTpath}/${uniq_id}_other.all_candsites -R ${ref_fa} -sr ${simple_repeat} > ${OUTpath}/${uniq_id}.other_sites.log 2>&1;" >> ${OUTPATH_ROOT}/05_cand_sites/${UNIQ_ID}_05candSites.sh
# done

# end
echo -e "job_begin\n name ${UNIQ_ID}_END\n memory 3G\n directory ${OUTPATH_ROOT}\n cmd_begin" >> ${sjmname}
echo -e "    /usr/bin/date;" >> ${sjmname}
#echo -e "    rm ${OUTPATH_ROOT}/01_porechop/*.trim.fastq.gz;" >> ${sjmname}
#echo -e "    rm ${OUTPATH_ROOT}/02_Nanofilt/*.trim.nanofilt.fastq.gz;" >> ${sjmname}
echo -e " cmd_end\njob_end\n" >> ${sjmname}

# total analysis order
# for FASTQ in ${INPATH}/*.fastq;
# do
#     uniq_id=${UNIQ_ID}_`basename ${FASTQ%.*}`
#     echo -e "order ${uniq_id}_00Nanoplot before ${UNIQ_ID}_00Nanoplot" >> ${sjmname}
# done
# for FASTQ in ${INPATH}/*.fastq;
# do
#     uniq_id=${UNIQ_ID}_`basename ${FASTQ%.*}`
#     echo -e "order ${UNIQ_ID}_00Nanoplot before ${uniq_id}_01porechop" >> ${sjmname}
#     echo -e "order ${uniq_id}_01porechop before ${UNIQ_ID}_01nanoplot" >> ${sjmname}
# done
# for FASTQ in ${INPATH}/*.fastq;
# do
#     uniq_id=${UNIQ_ID}_`basename ${FASTQ%.*}`
#     echo -e "order ${UNIQ_ID}_01nanoplot before ${uniq_id}_02Nanofilt" >> ${sjmname}
#     echo -e "order ${uniq_id}_02Nanofilt before ${UNIQ_ID}_02nanoplot" >> ${sjmname}
# done
# for FASTQ in ${INPATH}/*.fastq;
# do
#     uniq_id=${UNIQ_ID}_`basename ${FASTQ%.*}`
#     echo -e "order ${UNIQ_ID}_02nanoplot before ${uniq_id}_03minimap2" >> ${sjmname}
#     echo -e "order ${uniq_id}_03minimap2 before ${UNIQ_ID}_03nanoplot" >> ${sjmname}
# done
# for FASTQ in ${INPATH}/*.fastq;
# do
#     uniq_id=${UNIQ_ID}_`basename ${FASTQ%.*}`
#     echo -e "order ${UNIQ_ID}_03nanoplot before ${uniq_id}_04region" >> ${sjmname}
#     echo -e "order ${uniq_id}_04region before ${uniq_id}_05candSites" >> ${sjmname}
#     echo -e "order ${uniq_id}_05candSites before ${UNIQ_ID}_END" >> ${sjmname}
# done

for FASTQ in ${INPATH}/*.${SUFFIX_ORI}*;
do
    uniq_id=${UNIQ_ID}_`basename ${FASTQ} | sed "s/.${SUFFIX_ORI}//g"`
    echo -e "order ${UNIQ_ID}_START before ${uniq_id}_00Nanoplot" >> ${sjmname}
    # echo -e "order ${uniq_id}_00Nanoplot before ${uniq_id}_01porechop" >> ${sjmname}
    # echo -e "order ${uniq_id}_01porechop before ${uniq_id}_02Nanofilt" >> ${sjmname}
    echo -e "order ${uniq_id}_00Nanoplot before ${uniq_id}_02Nanofilt" >> ${sjmname}
    echo -e "order ${uniq_id}_02Nanofilt before ${uniq_id}_03minimap2" >> ${sjmname}
    if [[ ${LEVEL} == '1' ]]; then
        echo -e "order ${uniq_id}_03minimap2 before ${UNIQ_ID}_END" >> ${sjmname}
    else
        for key in ${chr_list[@]};
        do
            echo -e "order ${uniq_id}_03minimap2 before ${uniq_id}_04_bcftools_call" >> ${sjmname}
            echo -e "order ${uniq_id}_04_bcftools_call before ${uniq_id}_chr${key}_04_mpileup_call" >> ${sjmname}
            
            echo -e "order ${uniq_id}_chr${key}_04_mpileup_call before ${uniq_id}_04_callSite_CAT" >> ${sjmname}
            echo -e "order ${uniq_id}_04_callSite_CAT before ${UNIQ_ID}_END" >> ${sjmname}
        done
    fi
done
# for FASTQ in ${INPATH}/*.fastq*;
# do
#     uniq_id=${UNIQ_ID}_`basename ${FASTQ%.*} | sed 's/.fastq//g'`
#     echo -e "order ${uniq_id}_01porechop before ${uniq_id}_02Nanofilt" >> ${sjmname}
# done
# for FASTQ in ${INPATH}/*.fastq*;
# do
#     uniq_id=${UNIQ_ID}_`basename ${FASTQ%.*} | sed 's/.fastq//g'`
#     echo -e "order ${uniq_id}_02Nanofilt before ${uniq_id}_03minimap2" >> ${sjmname}
# done
# for FASTQ in ${INPATH}/*.fastq*;
# do
#     uniq_id=${UNIQ_ID}_`basename ${FASTQ%.*} | sed 's/.fastq//g'`
#     if [[ ${LEVEL} == '1' ]]; then
#         echo -e "order ${uniq_id}_03minimap2 before ${UNIQ_ID}_END" >> ${sjmname}
#     else
#         for key in ${chr_list[@]};
#         do
#             echo -e "order ${uniq_id}_03minimap2 before ${uniq_id}_04_callSite_chr${key}" >> ${sjmname}
#             echo -e "order ${uniq_id}_04_callSite_chr${key} before ${uniq_id}_04_callSite_CAT" >> ${sjmname}
#             echo -e "order ${uniq_id}_04_callSite_CAT before ${UNIQ_ID}_END" >> ${sjmname}
#         done
#     fi
# done
# if [[ ${LEVEL} == '1' ]]; then
#     echo -e "order ${UNIQ_ID}_03minimap2_level before ${UNIQ_ID}_END" >> ${sjmname}
# fi
# for FASTQ in ${INPATH}/*.fastq;
# do
#     uniq_id=${UNIQ_ID}_`basename ${FASTQ%.*}`
#     echo -e "order ${uniq_id}_03minimap2 before ${uniq_id}_04region" >> ${sjmname}
#     echo -e "order ${uniq_id}_04region before ${uniq_id}_05candSites" >> ${sjmname}
#     echo -e "order ${uniq_id}_05candSites before ${UNIQ_ID}_END" >> ${sjmname}
# done
# for FASTQ in ${INPATH}/*.fastq;
# do
#     uniq_id=${UNIQ_ID}_`basename ${FASTQ%.*}`
#     echo -e "order ${UNIQ_ID}_03nanoplot before ${uniq_id}_04region" >> ${sjmname}
#     echo -e "order ${uniq_id}_04region before ${UNIQ_ID}_END" >> ${sjmname}
# done
echo -e "log_dir ${OUTPATH_ROOT}/logs" >> ${sjmname}
    
# sjm --max_running 10 ${sjmname}

