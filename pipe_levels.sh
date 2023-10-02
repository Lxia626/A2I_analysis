#! /usr/bin/env bash

INpath=$1 # PATH for ORIGINAL FASTQ
OUTpath=$2 # PATH for storing output files "/public/home/Songlab/xial/work/path/"
uniq_id=$3 # UNIQUE ID for the data
single=$4 # {'1': single end; "2": paired end}
# sjm=$4 # whether to run in SJM {"0" [not use SJM], "1" [use SJM]}, Otherwise, "0" by default

## annotate:
echo "NOTE that:"
echo "  All log information are store in corresponding *.out file for each process, no matter it's running using qsub or sh."

# 01_fastqc
if [[ -d ${OUTpath}/01_fastqc ]]; then rm -r ${OUTpath}/01_fastqc; fi
if [[ -f ${OUTpath}/01_fastqc.sh ]]; then rm ${OUTpath}/01_fastqc.sh; fi
if [[ -f ${OUTpath}/01_fastqc.out ]]; then rm ${OUTpath}/01_fastqc.out; fi
mkdir ${OUTpath}/01_fastqc
touch ${OUTpath}/01_fastqc.sh
echo -e "#! /usr/bin/bash\n" >> ${OUTpath}/01_fastqc.sh
echo -e "#$ -cwd" >> ${OUTpath}/01_fastqc.sh
echo -e "#$ -o ${OUTpath}/01_fastqc.out" >> ${OUTpath}/01_fastqc.sh
echo -e "#$ -pe smp 1" >> ${OUTpath}/01_fastqc.sh
echo -e "#$ -N ${uniq_id}_01fastqc" >> ${OUTpath}/01_fastqc.sh
echo -e "#$ -j y" >> ${OUTpath}/01_fastqc.sh
echo -e "#$ -sync y\n" >> ${OUTpath}/01_fastqc.sh
echo -e "/public/home/Songlab/songyl/softwares/download/FastQC/fastqc -t 10 -o ${OUTpath}/01_fastqc/ ${INpath}/*.fastq" >> ${OUTpath}/01_fastqc.sh
qsub ${OUTpath}/01_fastqc.sh

# 02_cutadapt
if [[ -d ${OUTpath}/02_cutadapt ]]; then rm -r ${OUTpath}/02_cutadapt; fi
if [[ -f ${OUTpath}/02_cutadapt.sh ]]; then rm ${OUTpath}/02_cutadapt.sh; fi
if [[ -f ${OUTpath}/02_cutadapt.out ]]; then rm ${OUTpath}/02_cutadapt.out; fi
mkdir ${OUTpath}/02_cutadapt
touch ${OUTpath}/02_cutadapt.sh 
echo -e "#!  /usr/bin/bash" >> ${OUTpath}/02_cutadapt.sh
echo -e "#$ -cwd" >> ${OUTpath}/02_cutadapt.sh
echo -e "#$ -o ${OUTpath}/02_cutadapt.out" >> ${OUTpath}/02_cutadapt.sh
echo -e "#$ -pe smp 1" >> ${OUTpath}/02_cutadapt.sh
echo -e "#$ -N ${uniq_id}_02cutadapt" >> ${OUTpath}/02_cutadapt.sh
echo -e "#$ -j y" >> ${OUTpath}/02_cutadapt.sh
echo -e "#$ -sync y\n" >> ${OUTpath}/02_cutadapt.sh
if [[ ${single} == '2' ]]; then
    echo -e "for FASTQ in ${INpath}/*_1.fastq; do" >> ${OUTpath}/02_cutadapt.sh
    echo -e "" >> ${OUTpath}/02_cutadapt.sh
    echo -e "    inf1=\${FASTQ}" >> ${OUTpath}/02_cutadapt.sh
    echo -e "    inf2=\${FASTQ%_*}_2.fastq" >> ${OUTpath}/02_cutadapt.sh
    echo -e "    outf1=02_cutadapt/\`basename \${inf1}\`" >> ${OUTpath}/02_cutadapt.sh
    echo -e "    outf2=02_cutadapt/\`basename \${inf1%_*}\`_2.fastq" >> ${OUTpath}/02_cutadapt.sh
    echo -e "    lg_cut=02_cutadapt/\`basename \${inf1%_*}\`.log" >> ${OUTpath}/02_cutadapt.sh
    echo -e "    touch \$lg_cut" >> ${OUTpath}/02_cutadapt.sh
    echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/cutadapt -j 50 -a AGATCGGAAGAGCAC -A AGATCGGAAGAGCAC -q 20,20 --trim-n --length-tag 'length=' -m 30 -e 0.25 -o \$outf1 -p \$outf2 \$inf1 \$inf2 > \$lg_cut" >> ${OUTpath}/02_cutadapt.sh
    echo -e "    " >> ${OUTpath}/02_cutadapt.sh
    echo -e "done" >> ${OUTpath}/02_cutadapt.sh
else
    echo -e "for FASTQ in ${INpath}/*.fastq; do" >> ${OUTpath}/02_cutadapt.sh
    echo -e "" >> ${OUTpath}/02_cutadapt.sh
    echo -e "    inf1=\${FASTQ}" >> ${OUTpath}/02_cutadapt.sh
    echo -e "    outf1=02_cutadapt/\`basename \${inf1}\`" >> ${OUTpath}/02_cutadapt.sh
    echo -e "    lg_cut=02_cutadapt/\`basename \${inf1%_*}\`.log" >> ${OUTpath}/02_cutadapt.sh
    echo -e "    touch \$lg_cut" >> ${OUTpath}/02_cutadapt.sh
    echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/cutadapt -j 50 -a AGATCGGAAGAGCAC -q 20,20 --trim-n --length-tag 'length=' -m 30 -e 0.25 -o \$outf1 \$inf1 > \$lg_cut" >> ${OUTpath}/02_cutadapt.sh
    echo -e "    " >> ${OUTpath}/02_cutadapt.sh
    echo -e "done" >> ${OUTpath}/02_cutadapt.sh
fi 
qsub ${OUTpath}/02_cutadapt.sh

# 03_fastqc
if [[ -d ${OUTpath}/03_fastqc ]]; then rm -r ${OUTpath}/03_fastqc; fi
if [[ -f ${OUTpath}/03_fastqc.sh ]]; then rm ${OUTpath}/03_fastqc.sh; fi
if [[ -f ${OUTpath}/03_fastqc.out ]]; then rm ${OUTpath}/03_fastqc.out; fi
mkdir ${OUTpath}/03_fastqc
touch ${OUTpath}/03_fastqc.sh
echo -e "#! /usr/bin/bash\n" >> ${OUTpath}/03_fastqc.sh
echo -e "#$ -cwd" >> ${OUTpath}/03_fastqc.sh
echo -e "#$ -o ${OUTpath}/03_fastqc.out" >> ${OUTpath}/03_fastqc.sh
echo -e "#$ -pe smp 1" >> ${OUTpath}/03_fastqc.sh
echo -e "#$ -N ${uniq_id}_03fastqc" >> ${OUTpath}/03_fastqc.sh
echo -e "#$ -j y" >> ${OUTpath}/03_fastqc.sh
echo -e "#$ -sync y\n" >> ${OUTpath}/03_fastqc.sh
echo -e "/public/home/Songlab/songyl/softwares/download/FastQC/fastqc -t 10 -o ${OUTpath}/03_fastqc/ ${OUTpath}/02_cutadapt/*.fastq" >> ${OUTpath}/03_fastqc.sh
qsub ${OUTpath}/03_fastqc.sh

# 04_STAR [can not run in qsub-way, or it will raise error with FATAL ERROR. So instead, use nohup-way]
if [[ -d ${OUTpath}/04_STAR/ ]]; then rm -r ${OUTpath}/04_STAR/; fi
if [[ -f ${OUTpath}/04_STAR.sh ]]; then rm ${OUTpath}/04_STAR.sh; fi
if [[ -f ${OUTpath}/04_STAR.out ]]; then rm $${OUTpath}/04_STAR.out; fi
mkdir ${OUTpath}/04_STAR/
touch ${OUTpath}/04_STAR.sh
echo -e "#! /usr/bin/bash\n" >> ${OUTpath}/04_STAR.sh
echo -e "#$ -cwd" >> ${OUTpath}/04_STAR.sh
echo -e "#$ -o ${OUTpath}/04_STAR.out" >> ${OUTpath}/04_STAR.sh
echo -e "#$ -pe smp 1" >> ${OUTpath}/04_STAR.sh
echo -e "#$ -N ${uniq_id}_04_STAR" >> ${OUTpath}/04_STAR.sh
echo -e "#$ -j y" >> ${OUTpath}/04_STAR.sh
echo -e "#$ -sync y\n" >> ${OUTpath}/04_STAR.sh
if [[ ${single} == '2' ]]; then
    echo -e "for FASTQ in ${OUTpath}/02_cutadapt/*_1.fastq;" >> ${OUTpath}/04_STAR.sh
    echo -e "do" >> ${OUTpath}/04_STAR.sh
    echo -e "    inf1=\${FASTQ}" >> ${OUTpath}/04_STAR.sh
    echo -e "    inf2=\${FASTQ%_*}_2.fastq" >> ${OUTpath}/04_STAR.sh
    echo -e "    outf=${OUTpath}/04_STAR/\`basename \${FASTQ%_*}\`\".\"" >> ${OUTpath}/04_STAR.sh
    echo -e "    lg=${OUTpath}/04_STAR/\`basename \${FASTQ%_*}\`\".log\"" >> ${OUTpath}/04_STAR.sh
    echo -e "    /public/usr/bin/STAR --runThreadN 40 --genomeDir /public/home/Songlab/songyl/database/genome/human/GencodeV19/STAR_index/with_GencodeV19_gtf/ --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 4 --readFilesIn \$inf1 \$inf2 --outFileNamePrefix \$outf > \$lg" >> ${OUTpath}/04_STAR.sh
    echo -e "done" >> ${OUTpath}/04_STAR.sh
else
    echo -e "for FASTQ in ${OUTpath}/02_cutadapt/*.fastq;" >> ${OUTpath}/04_STAR.sh
    echo -e "do" >> ${OUTpath}/04_STAR.sh
    echo -e "    inf1=\${FASTQ}" >> ${OUTpath}/04_STAR.sh
    echo -e "    outf=${OUTpath}/04_STAR/\`basename \${FASTQ%_*}\`\".\"" >> ${OUTpath}/04_STAR.sh
    echo -e "    lg=${OUTpath}/04_STAR/\`basename \${FASTQ%_*}\`\".log\"" >> ${OUTpath}/04_STAR.sh
    echo -e "    /public/usr/bin/STAR --runThreadN 40 --genomeDir /public/home/Songlab/songyl/database/genome/human/GencodeV19/STAR_index/with_GencodeV19_gtf/ --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 4 --readFilesIn \$inf1 --outFileNamePrefix \$outf > \$lg" >> ${OUTpath}/04_STAR.sh
    echo -e "done" >> ${OUTpath}/04_STAR.sh
fi
sh ${OUTpath}/04_STAR.sh > ${OUTpath}/04_STAR.out 2>&1

# 05_featureCounts
if [[ -d ${OUTpath}/05_featureCounts/ ]]; then rm -r ${OUTpath}/05_featureCounts/; fi
# 05_featureCounts/exons
if [[ -d ${OUTpath}/05_featureCounts/exons ]]; then rm -r ${OUTpath}/05_featureCounts/exons; fi
if [[ -f ${OUTpath}/05_featureCounts_exons.sh ]]; then rm ${OUTpath}/05_featureCounts_exons.sh; fi
if [[ -f ${OUTpath}/05_featureCounts_exons.out ]]; then rm ${OUTpath}/05_featureCounts_exons.out; fi
mkdir -p ${OUTpath}/05_featureCounts/exons
touch ${OUTpath}/05_featureCounts_exons.sh
echo -e "#! /usr/bin/bash\n" >> ${OUTpath}/05_featureCounts_exons.sh
echo -e "#$ -cwd" >> ${OUTpath}/05_featureCounts_exons.sh
echo -e "#$ -o ${OUTpath}/05_featureCounts_exons.out" >> ${OUTpath}/05_featureCounts_exons.sh
echo -e "#$ -pe smp 1" >> ${OUTpath}/05_featureCounts_exons.sh
echo -e "#$ -N ${uniq_id}_05_FCounts_E" >> ${OUTpath}/05_featureCounts_exons.sh
echo -e "#$ -j y" >> ${OUTpath}/05_featureCounts_exons.sh
echo -e "#$ -sync y\n" >> ${OUTpath}/05_featureCounts_exons.sh
echo -e "for BAM in ${OUTpath}/04_STAR/*.Aligned.sortedByCoord.out.bam; do" >> ${OUTpath}/05_featureCounts_exons.sh
echo -e "    inf=\${BAM}" >> ${OUTpath}/05_featureCounts_exons.sh
echo -e "    outf=${OUTpath}/05_featureCounts/exons/\`basename \${BAM} | cut -d. -f1\`\".txt\"" >> ${OUTpath}/05_featureCounts_exons.sh
echo -e "    lg=${OUTpath}/05_featureCounts/exons/\`basename \${BAM} | cut -d. -f1\`\".log\"" >> ${OUTpath}/05_featureCounts_exons.sh
echo -e "    /public/usr/bin/featureCounts -p -t exon -g exon_id -s 2 -a /public/home/Songlab/songyl/database/annotation/human/gencode.v19.annotation.gtf -o \$outf \$inf > \$lg" >> ${OUTpath}/05_featureCounts_exons.sh
echo -e "done" >> ${OUTpath}/05_featureCounts_exons.sh
qsub ${OUTpath}/05_featureCounts_exons.sh
# 05_featureCounts/genes
if [[ -d ${OUTpath}/05_featureCounts/genes ]]; then rm -r ${OUTpath}/05_featureCounts/genes; fi
if [[ -f ${OUTpath}/05_featureCounts_genes.sh ]]; then rm ${OUTpath}/05_featureCounts_genes.sh; fi
if [[ -f ${OUTpath}/05_featureCounts_genes.out ]]; then rm ${OUTpath}/05_featureCounts_genes.out; fi
mkdir -p ${OUTpath}/05_featureCounts/genes
touch ${OUTpath}/05_featureCounts_genes.sh
echo -e "#! /usr/bin/bash\n" >> ${OUTpath}/05_featureCounts_genes.sh
echo -e "#$ -cwd" >> ${OUTpath}/05_featureCounts_genes.sh
echo -e "#$ -o ${OUTpath}/05_featureCounts_genes.out" >> ${OUTpath}/05_featureCounts_genes.sh
echo -e "#$ -pe smp 1" >> ${OUTpath}/05_featureCounts_genes.sh
echo -e "#$ -N ${uniq_id}_05_FCounts_G" >> ${OUTpath}/05_featureCounts_genes.sh
echo -e "#$ -j y" >> ${OUTpath}/05_featureCounts_genes.sh
echo -e "#$ -sync y\n" >> ${OUTpath}/05_featureCounts_genes.sh
echo -e "for BAM in ${OUTpath}/04_STAR/*.Aligned.sortedByCoord.out.bam; do" >> ${OUTpath}/05_featureCounts_genes.sh
echo -e "    inf=\${BAM}" >> ${OUTpath}/05_featureCounts_genes.sh
echo -e "    outf=${OUTpath}/05_featureCounts/genes/\`basename \${BAM} | cut -d. -f1\`\".txt\"" >> ${OUTpath}/05_featureCounts_genes.sh
echo -e "    lg=${OUTpath}/05_featureCounts/genes/\`basename \${BAM} | cut -d. -f1\`\".log\"" >> ${OUTpath}/05_featureCounts_genes.sh
echo -e "    /public/usr/bin/featureCounts -p -t exon -g gene_id -s 2 -a /public/home/Songlab/songyl/database/annotation/human/gencode.v19.annotation.gtf -o \$outf \$inf > \$lg" >> ${OUTpath}/05_featureCounts_genes.sh
echo -e "done" >> ${OUTpath}/05_featureCounts_genes.sh
qsub ${OUTpath}/05_featureCounts_genes.sh

# # 06_callKnownSitesLevels
if [[ -d ${OUTpath}/06_callKnownSitesLevels ]]; then rm -r ${OUTpath}/06_callKnownSitesLevels; fi
if [[ -f ${OUTpath}/06_callKnownSitesLevels.sh ]]; then rm ${OUTpath}/06_callKnownSitesLevels.sh; fi
if [[ -f ${OUTpath}/06_callKnownSitesLevels.out ]]; then rm ${OUTpath}/06_callKnownSitesLevels.out; fi
mkdir ${OUTpath}/06_callKnownSitesLevels
touch ${OUTpath}/06_callKnownSitesLevels.sh
echo -e "#! /usr/bin/bash\n" >> ${OUTpath}/06_callKnownSitesLevels.sh
echo -e "#$ -cwd" >> ${OUTpath}/06_callKnownSitesLevels.sh
echo -e "#$ -o ${OUTpath}/06_callKnownSitesLevels.out" >> ${OUTpath}/06_callKnownSitesLevels.sh
echo -e "#$ -pe smp 1" >> ${OUTpath}/06_callKnownSitesLevels.sh
echo -e "#$ -N ${uniq_id}_06_KnownSitesLevels" >> ${OUTpath}/06_callKnownSitesLevels.sh
echo -e "#$ -j y" >> ${OUTpath}/06_callKnownSitesLevels.sh
echo -e "#$ -sync y\n" >> ${OUTpath}/06_callKnownSitesLevels.sh
echo -e "/usr/bin/perl /public/home/Songlab/xial/scripts/pipeline_songyl/call_levels/editing_site_call_from_RNAseq_ADAR_20170512.pl --in ${OUTpath}/04_STAR/ --out ${OUTpath}/06_callKnownSitesLevels" >> ${OUTpath}/06_callKnownSitesLevels.sh
qsub ${OUTpath}/06_callKnownSitesLevels.sh
# 06_countMappingRates
if [[ -d ${OUTpath}/06_countMappingRates ]]; then rm -r ${OUTpath}/06_countMappingRates; fi
if [[ -f ${OUTpath}/06_countMappingRates.sh ]]; then rm ${OUTpath}/06_countMappingRates.sh; fi
if [[ -f ${OUTpath}/06_countMappingRates.out ]]; then rm ${OUTpath}/06_countMappingRates.out; fi
mkdir ${OUTpath}/06_countMappingRates
touch ${OUTpath}/06_countMappingRates.sh
echo -e "#! /usr/bin/bash\n" >> ${OUTpath}/06_countMappingRates.sh
echo -e "#$ -cwd" >> ${OUTpath}/06_countMappingRates.sh
echo -e "#$ -o ${OUTpath}/06_countMappingRates.out" >> ${OUTpath}/06_countMappingRates.sh
echo -e "#$ -pe smp 1" >> ${OUTpath}/06_countMappingRates.sh
echo -e "#$ -N ${uniq_id}_06_MappingRates" >> ${OUTpath}/06_countMappingRates.sh
echo -e "#$ -j y" >> ${OUTpath}/06_countMappingRates.sh
echo -e "#$ -sync y\n" >> ${OUTpath}/06_countMappingRates.sh
echo -e "echo -e \"Sample\tUnique\tMultiple\tTotal\" > ${OUTpath}/06_countMappingRates/STAR_mapping_rates_Res" >> ${OUTpath}/06_countMappingRates.sh
echo -e "for LOG_OUT in ${OUTpath}/04_STAR/*.Log.final.out; do" >> ${OUTpath}/06_countMappingRates.sh
echo -e "    sample_id=\`basename \${LOG_OUT} | cut -d. -f1\`" >> ${OUTpath}/06_countMappingRates.sh
echo -e "    Uniq_ratio=\`grep 'Uniquely mapped reads %' \${LOG_OUT} | awk -F '\t' '{print \$2}' | cut -d% -f1\`" >> ${OUTpath}/06_countMappingRates.sh
echo -e "    Multi_ratio=\`grep '% of reads mapped to' \${LOG_OUT} | awk 'BEGIN{sum=0}{sum+=\$NF}END{print sum}'\`" >> ${OUTpath}/06_countMappingRates.sh
echo -e "    Total_ratio=\`echo \"scale=2;\${Multi_ratio}+\${Uniq_ratio}\" | bc\`" >> ${OUTpath}/06_countMappingRates.sh
echo -e "    echo -e \"\${sample_id}\t\${Uniq_ratio}%\t\${Multi_ratio}%\t\${Total_ratio}%\" >> ${OUTpath}/06_countMappingRates/STAR_mapping_rates_Res" >> ${OUTpath}/06_countMappingRates.sh
echo -e "done" >> ${OUTpath}/06_countMappingRates.sh
qsub ${OUTpath}/06_countMappingRates.sh

# 07_RNAEditingIndexer
if [[ -d ${OUTpath}/07_RNAEditingIndexer ]]; then rm -r ${OUTpath}/07_RNAEditingIndexer; fi
if [[ -f ${OUTpath}/07_RNAEditingIndexer.sh ]]; then rm ${OUTpath}/07_RNAEditingIndexer.sh; fi
if [[ -f ${OUTpath}/07_RNAEditingIndexer.out ]]; then rm ${OUTpath}/07_RNAEditingIndexer.out; fi
mkdir ${OUTpath}/07_RNAEditingIndexer
touch ${OUTpath}/07_RNAEditingIndexer.sh
echo -e "#! /usr/bin/bash\n" >> ${OUTpath}/07_RNAEditingIndexer.sh
echo -e "#$ -cwd" >> ${OUTpath}/07_RNAEditingIndexer.sh
echo -e "#$ -o ${OUTpath}/07_RNAEditingIndexer.out" >> ${OUTpath}/07_RNAEditingIndexer.sh
echo -e "#$ -pe smp 1" >> ${OUTpath}/07_RNAEditingIndexer.sh
echo -e "#$ -N ${uniq_id}_07_RNAEditingIndexer" >> ${OUTpath}/07_RNAEditingIndexer.sh
echo -e "#$ -j y" >> ${OUTpath}/07_RNAEditingIndexer.sh
echo -e "#$ -sync y\n" >> ${OUTpath}/07_RNAEditingIndexer.sh
echo -e "/public/home/Songlab/songyl/softwares/RNAEditingIndexer-master/src/RNAEditingIndex/RNAEditingIndex -d ${OUTpath}/04_STAR/ -f Aligned.sortedByCoord.out.bam -l ${OUTpath}/07_RNAEditingIndexer/default_log/ -o ${OUTpath}/07_RNAEditingIndexer/default_o/ -os ${OUTpath}/07_RNAEditingIndexer/default_os/ --genome hg19 --paired_end" >> ${OUTpath}/07_RNAEditingIndexer.sh
sh ${OUTpath}/07_RNAEditingIndexer.sh > ${OUTpath}/07_RNAEditingIndexer.out 2>&1

#### running
#qsub -o ${OUTpath}/01_fastqc.out -pe smp 1 -j y -cwd -b y -N ${uniq_id}_01 -sync y ${OUTpath}/01_fastqc.sh
#qsub -o ${OUTpath}/02_cutadapt.out -j y -pe smp 1 -cwd -b y -N ${uniq_id}_02 -sync y ${OUTpath}/02_cutadapt.sh
#qsub -o ${OUTpath}/03_fastqc.out -j y -pe smp 1 -cwd -b y -N ${uniq_id}_03 -sync y ${OUTpath}/03_fastqc.sh
#sh ${OUTpath}/04_STAR.sh > ${OUTpath}/04_STAR.out 2>&1
#qsub -o ${OUTpath}/05_featureCounts_exons.out -pe smp 1 -j y -cwd -b y -N ${uniq_id}_05exons -sync y ${OUTpath}/05_featureCounts_exons.sh
#qsub -o ${OUTpath}/05_featureCounts_genes.out -pe smp 1 -j y -cwd -b y -N ${uniq_id}_05genes -sync y ${OUTpath}/05_featureCounts_genes.sh
#qsub -o ${OUTpath}/06_callKnownSitesLevels.out -pe smp 1 -j y -cwd -b y -N ${uniq_id}_06levels -sync y ${OUTpath}/06_callKnownSitesLevels.sh
#qsub -o ${OUTpath}/06_countMappingRates.out -pe smp 1 -j y -cwd -b y -N ${uniq_id}_06rates -sync y ${OUTpath}/06_countMappingRates.sh
#sh ${OUTpath}/07_RNAEditingIndexer.sh > ${OUTpath}/07_RNAEditingIndexer.out 2>&1


# #### pipeline in SJM
# if [[ $sjm -eq 1 ]]; then

#     SJMFile=${OUTpath}/${uniq_id}".levels.sjm"
#     if [[ -f ${SJMFile} ]]; then rm $SJMFile; fi
#     touch $SJMFile
#     echo -e "job_begin\n name ${uniq_id}_01_fastqc\n memory 10G\n directory ${OUTpath}\n cmd_begin" >> $SJMFile
#     echo -e "/usr/bin/bash ${OUTpath}/01_fastqc.sh > ${OUTpath}/01_fastqc.log 2>&1;" >> $SJMFile
#     echo -e "cmd_end\njob_end\n" >> $SJMFile

#     echo -e "job_begin\n name ${uniq_id}_02_cutadapt\n memory 30G\n directory ${OUTpath}\n cmd_begin" >> $SJMFile
#     echo -e "/usr/bin/bash ${OUTpath}/02_cutadapt.sh > ${OUTpath}/02_cutadapt.log 2>&1;" >> $SJMFile
#     echo -e "cmd_end\njob_end\n" >> $SJMFile

#     echo -e "job_begin\n name ${uniq_id}_03_fastqc\n memory 10G\n directory ${OUTpath}\n cmd_begin" >> $SJMFile
#     echo -e "/usr/bin/bash ${OUTpath}/03_fastqc.sh > ${OUTpath}/03_fastqc.log 2>&1;" >> $SJMFile
#     echo -e "cmd_end\njob_end\n" >> $SJMFile

#     echo -e "job_begin\n name ${uniq_id}_04_STAR\n memory 50G\n directory ${OUTpath}\n cmd_begin" >> $SJMFile
#     echo -e "/usr/bin/bash ${OUTpath}/04_STAR.sh > ${OUTpath}/04_STAR.log 2>&1;" >> $SJMFile
#     echo -e "cmd_end\njob_end\n" >> $SJMFile

#     echo -e "job_begin\n name ${uniq_id}_05_featureCounts\n memory 30G\n directory ${OUTpath}\n cmd_begin" >> $SJMFile
#     echo -e "/usr/bin/bash ${OUTpath}/05_featureCounts_exons.sh > ${OUTpath}/05_featureCounts_exons.log 2>&1;" >> $SJMFile
#     echo -e "/usr/bin/bash ${OUTpath}/05_featureCounts_genes.sh > ${OUTpath}/05_featureCounts_genes.log 2>&1;" >> $SJMFile
#     echo -e "cmd_end\njob_end\n" >> $SJMFile

#     echo -e "job_begin\n name ${uniq_id}_06_knownSites\n memory 10G\n directory ${OUTpath}\n cmd_begin" >> $SJMFile
#     echo -e "/usr/bin/bash ${OUTpath}/06_callKnownSitesLevels.sh > ${OUTpath}/06_callKnownSitesLevels.log 2>&1;" >> $SJMFile
#     echo -e "cmd_end\njob_end\n" >> $SJMFile

#     echo -e "job_begin\n name ${uniq_id}_06_mappingRates\n memory 10G\n directory ${OUTpath}\n cmd_begin" >> $SJMFile
#     echo -e "/usr/bin/bash ${OUTpath}/06_countMappingRates.sh > ${OUTpath}/06_mappingRates.log 2>&1;" >> $SJMFile
#     echo -e "cmd_end\njob_end\n" >> $SJMFile

#     echo -e "job_begin\n name ${uniq_id}_07_RNAEditingIndexer\n memory 30G\n directory ${OUTpath}\n cmd_begin" >> $SJMFile
#     echo -e "/usr/bin/bash ${OUTpath}/07_RNAEditingIndexer.sh > ${OUTpath}/07_RNAEditingIndexer.log 2>&1;" >> $SJMFile
#     echo -e "cmd_end\njob_end\n" >> $SJMFile

#     if [[ -d ${OUTpath}/logs ]]; then rm -r ${OUTpath}/logs; fi
#     mkdir ${OUTpath}/logs
#     echo -e "order ${uniq_id}_01_fastqc before ${uniq_id}_02_cutadapt" >> ${SJMFile}
#     echo -e "order ${uniq_id}_02_cutadapt before ${uniq_id}_03_fastqc" >> ${SJMFile}
#     echo -e "order ${uniq_id}_03_fastqc before ${uniq_id}_04_STAR" >> ${SJMFile}
#     echo -e "order ${uniq_id}_04_STAR before ${uniq_id}_05_featureCounts" >> ${SJMFile}
#     echo -e "order ${uniq_id}_05_featureCounts before ${uniq_id}_06_knownSites" >> ${SJMFile}
#     echo -e "order ${uniq_id}_05_featureCounts before ${uniq_id}_06_mappingRates" >> ${SJMFile}
#     echo -e "order ${uniq_id}_06_knownSites before ${uniq_id}_07_RNAEditingIndexer" >> ${SJMFile}
#     echo -e "order ${uniq_id}_06_mappingRates before ${uniq_id}_07_RNAEditingIndexer" >> ${SJMFile}
#     echo -e "log_dir ${OUTpath}/logs" >> ${SJMFile}

# fi
