#! /usr/bin/env bash

# not need the absolute path
INPATH=`pwd`/$1
OUTPATH_ROOT=`pwd`/$2
#SEQUENCING_SUMMARY=$3 # the sequencing summary file from albacore
#FAST5_FILE=$4 # the raw ONT signal files
UNIQ_ID=$3
LENGTH=$4 # if LENGTH==0, then all reads were maintained
QUAL=$5 # limit of base quality
DEPTH=$6 # minimum depth of all candidate sites
DELETE_ORIGINAL=$7 #{'0', '1'}

# basic settings
# GRCh37_75_genome="/public/home/Songlab/songyl/database/genome/human/GRCh37_release75/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
REF_FA=/public/home/Songlab/xial/Work/database/hg19/hg19_softmasked.fa
GTF_FILE=/public/home/Songlab/xial//Work/database/hg19/gencode.v27lift37.annotation.sorted.gtf.gz
ALU_BED=/public/home/Songlab/songyl/database/hg19_Alu.bed
NONALU_BED=/public/home/Songlab/songyl/database/Repeat/human/human_nonAlu.bed

if [[ ! -d ${OUTPATH_ROOT} ]]; then 
    mkdir ${OUTPATH_ROOT}
fi
# 00: post-run basecalling with the lastest version of Guppy after live basecalling using an older version of Guppy in MinKNOW.
# The process was based on Zhong et al., Nat Commun, 2023
if [[ ${DELETE_ORIGINAL} == '1' ]]; then
    if [[ -d ${OUTPATH_ROOT}/00_Nanoplot ]]; then rm -r ${OUTPATH_ROOT}/00_Nanoplot; fi
    mkdir -p ${OUTPATH_ROOT}/00_Nanoplot
fi
if [[ -f ${OUTPATH_ROOT}/00_Nanoplot/${UNIQ_ID}_00Nanoplot.sh ]]; then rm ${OUTPATH_ROOT}/00_Nanoplot/${UNIQ_ID}_00Nanoplot.sh; fi
touch ${OUTPATH_ROOT}/00_Nanoplot/${UNIQ_ID}_00Nanoplot.sh
for FASTQ in ${INPATH}/*.fastq;
do
    uniq_id=${UNIQ_ID}_`basename ${FASTQ%.*}`
    echo -e "echo \"Analyzing ${uniq_id}_00Nanoplot\"" >> ${OUTPATH_ROOT}/00_Nanoplot/${UNIQ_ID}_00Nanoplot.sh
    echo -e "/public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoPlot -t 8 --fastq ${FASTQ} --loglength -o ${OUTPATH_ROOT}/00_Nanoplot --format pdf -p ${uniq_id}. --plots dot hex > ${OUTPATH_ROOT}/00_Nanoplot/${uniq_id}_nanoplot.log 2>&1 &" >> ${OUTPATH_ROOT}/00_Nanoplot/${UNIQ_ID}_00Nanoplot.sh
done
echo -e "echo \"Analyzing ${UNIQ_ID}_00Nanoplot\"" >> ${OUTPATH_ROOT}/00_Nanoplot/${UNIQ_ID}_00Nanoplot.sh
# echo -e "/public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoPlot -t 8 --fastq ${INPATH}/*.fastq --loglength -o ${OUTPATH_ROOT}/00_Nanoplot --format pdf -p ${UNIQ_ID}_nanoplot. --plots dot hex > ${OUTPATH_ROOT}/00_Nanoplot/nanoplot.log 2>&1 &" >> ${OUTPATH_ROOT}/00_Nanoplot/${UNIQ_ID}_00Nanoplot.sh
# echo -e "wait" >> ${OUTPATH_ROOT}/00_Nanoplot/${UNIQ_ID}_00Nanoplot.sh

# 01: porechop
if [[ ${DELETE_ORIGINAL} == '1' ]]; then
    if [[ -d ${OUTPATH_ROOT}/01_porechop ]]; then rm -r ${OUTPATH_ROOT}/01_porechop; fi
    mkdir -p ${OUTPATH_ROOT}/01_porechop
fi
if [[ -f ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh ]]; then rm ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh; fi
touch ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh
echo -e "porechop_process() {" >> ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh
echo -e "    FASTQ=\$1" >> ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh
echo -e "    OUTpath=\$2" >> ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh
echo -e "    uniq_id=\$3" >> ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh
echo -e "    echo \"Analyzing ${uniq_id}_01_porechop\"" >> ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh
echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/porechop -i \${FASTQ} -o \${OUTpath}/\${uniq_id}.trim.fastq.gz --discard_middle --format fastq.gz > \${OUTpath}/\${uniq_id}.porechop.log 2>&1" >> ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh
echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoPlot -t 8 --fastq \${OUTpath}/\${uniq_id}.trim.fastq.gz --loglength -o \${OUTpath}/nanoplot --format pdf -p \${uniq_id}. --plots dot hex > \${OUTpath}/\${uniq_id}_nanoplot.log 2>&1" >> ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh
echo -e "}" >> ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh
echo -e "\n\n" >> ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh
for FASTQ in ${INPATH}/*.fastq;
do
    uniq_id=${UNIQ_ID}_`basename ${FASTQ%.*}`
    OUTpath=${OUTPATH_ROOT}/01_porechop
    echo -e "porechop_process ${FASTQ} ${OUTpath} ${uniq_id} &" >> ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh
done
echo -e "wait" >> ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh
# echo -e "echo \"Analyzing ${uniq_id}_01_nanoplot\"" >> ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh
# echo -e "/public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoPlot -t 8 --fastq ${OUTpath}/${UNIQ_ID}*.trim.fastq.gz --loglength -o ${OUTpath}/nanoplot --format pdf -p ${UNIQ_ID}_nanoplot. --plots dot hex > ${OUTpath}/nanoplot.log 2>&1" >> ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh

# 02: Nanofilt
if [[ ${DELETE_ORIGINAL} == '1' ]]; then
    if [[ -d ${OUTPATH_ROOT}/02_Nanofilt ]]; then rm -r ${OUTPATH_ROOT}/02_Nanofilt; fi
    mkdir -p ${OUTPATH_ROOT}/02_Nanofilt
fi
if [[ -f ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh ]]; then rm ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh; fi
touch ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
echo -e "nanofilt_process() {" >> ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
echo -e "    OUTPATH_ROOT=\$1" >> ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
echo -e "    uniq_id=\$2" >> ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
echo -e "    OUTpath=\$3" >> ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
echo -e "    echo \"Analyzing ${uniq_id}_02_Nanofilt\"" >> ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
if [[ ${LENGTH} -eq 0 ]]; then
    echo -e "    gunzip -c \${OUTPATH_ROOT}/01_porechop/\${uniq_id}.trim.fastq.gz | /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoFilt -q ${QUAL} --headcrop 50 | gzip > \${OUTpath}/\${uniq_id}.trim.nanofit.fastq.gz" >> ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
else
    echo -e "    gunzip -c \${OUTPATH_ROOT}/01_porechop/\${uniq_id}.trim.fastq.gz | /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoFilt -l ${LENGTH} -q ${QUAL} --headcrop 50 | gzip > \${OUTpath}/\${uniq_id}.trim.nanofit.fastq.gz" >> ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
fi
echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoPlot -t 8 --fastq \${OUTpath}/\${uniq_id}.trim.nanofit.fastq.gz --loglength -o \${OUTpath}/nanoplot --format pdf -p \${uniq_id}. --plots dot hex > \${OUTpath}/\${uniq_id}.nanoplot.log 2>&1" >> ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
echo -e "}" >> ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
echo -e "\n\n" >> ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
for FASTQ in ${INPATH}/*.fastq;
do
    uniq_id=${UNIQ_ID}_`basename ${FASTQ%.*}`
    OUTpath=${OUTPATH_ROOT}/02_Nanofilt
    echo -e "nanofilt_process ${OUTPATH_ROOT} ${uniq_id} ${OUTpath} &" >> ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh 
    done
echo -e "wait" >> ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
# echo -e "echo \"Analyzing ${uniq_id}_02_nanoplot\"" >> ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
# echo -e "/public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoPlot -t 8 --fastq ${OUTpath}/${UNIQ_ID}*.trim.nanofit.fastq.gz --loglength -o ${OUTpath}/nanoplot --format pdf -p ${UNIQ_ID}_nanoplot. --plots dot hex > ${OUTpath}/nanoplot.log 2>&1" >> ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh
    
# 03: minimap2 and editing sites calling
if [[ ${DELETE_ORIGINAL} == '1' ]]; then
    if [[ -d ${OUTPATH_ROOT}/03_minimap2/ ]]; then rm -r ${OUTPATH_ROOT}/03_minimap2/; fi
    mkdir -p ${OUTPATH_ROOT}/03_minimap2/
fi
if [[ -f ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03minimap2.sh ]]; then rm ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03minimap2.sh; fi
touch ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03minimap2.sh
echo -e "minimap2_process() {" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03minimap2.sh
echo -e "    OUTPATH_ROOT=\$1" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03minimap2.sh
echo -e "    uniq_id=\$2" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03minimap2.sh
echo -e "    OUTpath=\$3" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03minimap2.sh
echo -e "    echo \"Analyzing ${uniq_id}_03_minimap2\"" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03minimap2.sh
echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/minimap2 -t 10 --secondary=no --cs -ax splice -uf -k14 ${REF_FA} \${OUTPATH_ROOT}/02_Nanofilt/\${uniq_id}.trim.nanofit.fastq.gz | /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/samtools view -q 30 -b -@ 4 -F 2052 | /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/samtools sort -@ 4 -O BAM -o \${OUTpath}/\${uniq_id}.bam > \${OUTpath}/\${uniq_id}.log 2>&1"  >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03minimap2.sh # parameter '--cs' and '-F 2052' were followed by Liu et al., Genome Biol, 2023; parameter '-ax splice -uf -k14' was followed by Workman et al., Nat Methods, 2019
echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/samtools index \${OUTpath}/\${uniq_id}.bam" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03minimap2.sh
echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoPlot -t 8 --bam \${OUTpath}/\${uniq_id}.bam --loglength -o \${OUTpath}/nanoplot --format pdf -p \${uniq_id}. --plots dot hex > \${OUTpath}/\${uniq_id}.nanoplot.log 2>&1" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03minimap2.sh
echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/samtools flagstat \${OUTpath}/\${uniq_id}.bam > \${OUTpath}/\${uniq_id}.bam.flagstat" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03minimap2.sh
echo -e "    /public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/python /public/home/Songlab/xial/Work/software/L-GIREMI/bin/correct_read_strand -b \{OUTpath}/\${uniq_id}.bam -o \${OUTpath}/\${uniq_id}.bam -t 10 --annotation_gtf ${GTF_FILE} --genome_fasta ${REF_FA} > \${OUTpath}/\${uniq_id}.correctStrand.log 2>&1" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03minimap2.sh
echo -e "    /public/usr/bin/bcftools mpileup -B -q ${QUAL} -Q 0 -I -d 10000 --threads 10 -f ${REF_FA} \${OUTpath}/\${uniq_id}.bam | /public/usr/bin/bcftools call -P0.01 --threads 10 --ploidy 1 -vm -Oz -o \${OUTpath}/\${uniq_id}.vcf.gz > \${OUTpath}/\${uniq_id}.bcftools.log 2>&1" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03minimap2.sh # -Q7 may impacts the DP4 parameter. Parameter '-I', '-Q', '--ploidy' and '-m' were used according to Hall et al., Lancet Microbe, 2022. '-I' may hide some homo-edited sites?
echo -e "    /public/usr/bin/bcftools index \${OUTpath}/\${uniq_id}.vcf.gz" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03minimap2.sh
echo -e "}" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03minimap2.sh
echo -e "\n\n" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03minimap2.sh
for FASTQ in ${INPATH}/*.fastq;
do
    uniq_id=${UNIQ_ID}_`basename ${FASTQ%.*}`
    OUTpath=${OUTPATH_ROOT}/03_minimap2
    echo -e "minimap2_process ${OUTPATH_ROOT} ${uniq_id} ${OUTpath} &" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03minimap2.sh
done
echo -e "wait" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03minimap2.sh
# echo -e "echo \"Analyzing ${uniq_id}_03_nanoplot\"" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03minimap2.sh
# echo -e "/public/home/Songlab/xial/miniconda3/envs/Nanopore/bin/NanoPlot -t 8 --bam ${OUTpath}/${UNIQ_ID}*.bam --loglength -o ${OUTpath}/nanoplot --format pdf -p ${UNIQ_ID}_nanoplot. --plots dot hex > ${OUTpath}/nanoplot.log 2>&1" >> ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03minimap2.sh
    
# 04: dividing sites into Alu, nonAlu-repeat, other
snp_db="/public/home/Songlab/xial/Work/database/hg19/SNP_all_dbSNP_1000Genomes_UWash.txt"
if [[ ${DELETE_ORIGINAL} == '1' ]]; then
    if [[ -d ${OUTPATH_ROOT}/04_region/ ]]; then rm -r ${OUTPATH_ROOT}/04_region/; fi
    mkdir -p ${OUTPATH_ROOT}/04_region/
fi
count_dbSNP_depth() {
    IN_VCF=$1
    sjmname=$2
    OUTpath=$3
    uniq_id=$4
    region=$5
    echo -e "    original=\`less ${IN_VCF} | grep -v '##' | wc -l\`;" >> ${sjmname}
    echo -e "    /public/usr/bin/bcftools filter --threads 4 --mask-file ~/Work/database/hg19/SNP_all_dbSNP_1000Genomes_UWash.txt -s SNP ${IN_VCF} | grep -v 'SNP' | gzip -c > ${OUTpath}/${uniq_id}.${region}.dbSNP_filt.vcf.gz;" >> ${sjmname}
    echo -e "    dbSNP_filt=\`less ${OUTpath}/${uniq_id}.${region}.dbSNP_filt.vcf.gz | grep -v '##' | wc -l\`;" >> ${sjmname}
    echo -e "    /public/usr/bin/bcftools filter --threads 4 --exclude 'DP4[0]+DP4[1]+DP4[2]+DP4[3]<10' -s DEPTH_FAILED ${OUTpath}/${uniq_id}.${region}.dbSNP_filt.vcf.gz | grep -v 'DEPTH_FAILED' | gzip -c > ${OUTpath}/${uniq_id}.${region}.dbSNP_depth_filt.vcf.gz;" >> ${sjmname}
    echo -e "    depth_filt=\`less ${OUTpath}/${uniq_id}.${region}.dbSNP_depth_filt.vcf.gz | grep -v '##' | wc -l\`;" >> ${sjmname}
    echo -e "    echo -e \"dbSNP number for `basename ${IN_VCF}`: \`expr \${original} - \${dbSNP_filt}\`\" >> ${OUTpath}/${uniq_id}.count;" >> ${sjmname}
    echo -e "    echo -e \"depth filtered number for `basename ${IN_VCF}`: \`expr \${dbSNP_filt} - \${depth_filt}\`\" >> ${OUTpath}/${uniq_id}.count;" >> ${sjmname}
}
if [[ -f ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh ]]; then rm ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh; fi
touch ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
echo -e "region_process() {" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
echo -e "    OUTPATH_ROOT=\$1" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
echo -e "    uniq_id=\$2" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
echo -e "    OUTpath=\$3" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
echo -e "    prefix=\${OUTpath}/\${uniq_id}" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
echo -e "    echo \"Analyzing ${uniq_id}_04_region\"" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
# exclude SNPs
echo -e "    /public/home/Songlab/xial/miniconda3/bin/vcftools --exclude-positions ${snp_db} --gzvcf \${OUTPATH_ROOT}/03_minimap2/\${uniq_id}.vcf.gz --recode --recode-INFO-all --stdout | /public/usr/bin/bcftools view -Oz -o \${prefix}.snpdb.vcf.gz > \${prefix}.ALU.snpdb.log 2>&1" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
# ALU region
echo -e "    /public/home/Songlab/xial/miniconda3/bin/vcftools --bed ${ALU_BED} --gzvcf \${prefix}.snpdb.vcf.gz --recode-INFO-all --recode --stdout | /public/usr/bin/bcftools view -Oz -o \${prefix}.ALU.vcf.gz > \${prefix}.ALU.gz.log 2>&1" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
echo -e "    /public/usr/bin/bcftools filter --threads 4 --include 'DP4[0]+DP4[1]+DP4[2]+DP4[3]>${DEPTH}' -o \${prefix}.ALU.filtered.bcf -Ob \${prefix}.ALU.vcf.gz > \${prefix}.ALU.filtered.log 2>&1" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
echo -e "    /public/usr/bin/bcftools index \${prefix}.ALU.filtered.bcf" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh 
echo -e "    /public/home/Songlab/xial/miniconda3/bin/vcftools --exclude-bed ${ALU_BED} --gzvcf \${prefix}.snpdb.vcf.gz --recode-INFO-all --recode --out \${prefix}.nonALU > \${prefix}.nonALU.log 2>&1" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
# nonALU repeat region
echo -e "    /public/home/Songlab/xial/miniconda3/bin/vcftools --bed ${NONALU_BED} --vcf \${prefix}.nonALU.recode.vcf --recode-INFO-all --recode --stdout | /public/usr/bin/bcftools view -Oz -o \${prefix}.nonALU_rep.vcf.gz > \${prefix}.nonALU_rep.gz.log 2>&1" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
echo -e "    /public/usr/bin/bcftools filter --threads 4 --include 'DP4[0]+DP4[1]+DP4[2]+DP4[3]>${DEPTH}' -o \${prefix}.nonALU_rep.filtered.bcf -Ob \${prefix}.nonALU_rep.vcf.gz > \${prefix}.nonALU_rep.filtered.log 2>&1" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
echo -e "    /public/usr/bin/bcftools index \${prefix}.nonALU_rep.filtered.bcf" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
# other
echo -e "    /public/home/Songlab/xial/miniconda3/bin/vcftools --exclude-bed ${NONALU_BED} --vcf \${prefix}.nonALU.recode.vcf --recode-INFO-all --recode --stdout | /public/usr/bin/bcftools view -Oz -o \${prefix}.other.vcf.gz > \${prefix}.other.gz.log 2>&1" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
echo -e "    /public/usr/bin/bcftools filter --threads 4 --include 'DP4[0]+DP4[1]+DP4[2]+DP4[3]>${DEPTH}' -o \${prefix}.other.filtered.bcf -Ob \${prefix}.other.vcf.gz > \${prefix}.other.filtered.log 2>&1" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
echo -e "    /public/usr/bin/bcftools index \${prefix}.other.filtered.bcf" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
# calculate
echo -e "    ALU_num=\`/public/usr/bin/bcftools view \${prefix}.ALU.filtered.bcf | grep -v '##' | wc -l\`" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
echo -e "    nonALU_rep_num=\`/public/usr/bin/bcftools view \${prefix}.nonALU_rep.filtered.bcf | grep -v '##' | wc -l\`" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
echo -e "    other_num=\`/public/usr/bin/bcftools view \${prefix}.other.filtered.bcf | grep -v '##' | wc -l\`" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
echo -e "    echo -e \"ALU region number: \${ALU_num}\" >> \${prefix}.count" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
echo -e "    echo -e \"nonALU_rep region number: \${nonALU_rep_num}\" >> \${prefix}.count" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
echo -e "    echo -e \"Other region number: \${other_num}\" >> \${prefix}.count" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
echo -e "}" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
echo -e "\n\n" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
for FASTQ in ${INPATH}/*.fastq;
do
    uniq_id=${UNIQ_ID}_`basename ${FASTQ%.*}`
    OUTpath=${OUTPATH_ROOT}/04_region
    echo -e "region_process ${OUTPATH_ROOT} ${uniq_id} ${OUTpath} &" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh
done
echo -e "wait" >> ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh

# 05 candidate site filtration
simple_repeat="/public/home/Songlab/songyl/database/hg19_SimpleRepeat.txt"
if [[ ${DELETE_ORIGINAL} == '1' ]]; then
    if [[ -d ${OUTPATH_ROOT}/05_cand_sites/ ]]; then rm -r ${OUTPATH_ROOT}/05_cand_sites/; fi
    mkdir -p ${OUTPATH_ROOT}/05_cand_sites/
fi
if [[ -f ${OUTPATH_ROOT}/05_cand_sites/${UNIQ_ID}_05candSites.sh ]]; then rm ${OUTPATH_ROOT}/05_cand_sites/${UNIQ_ID}_05candSites.sh; fi
touch ${OUTPATH_ROOT}/05_cand_sites/${UNIQ_ID}_05candSites.sh
echo -e "candSites_process() {" >> ${OUTPATH_ROOT}/05_cand_sites/${UNIQ_ID}_05candSites.sh
echo -e "    OUTPATH_ROOT=\$1" >> ${OUTPATH_ROOT}/05_cand_sites/${UNIQ_ID}_05candSites.sh
echo -e "    uniq_id=\$2" >> ${OUTPATH_ROOT}/05_cand_sites/${UNIQ_ID}_05candSites.sh
echo -e "    OUTpath=\$3" >> ${OUTPATH_ROOT}/05_cand_sites/${UNIQ_ID}_05candSites.sh
echo -e "    echo \"Analyzing ${uniq_id}_05_cand_sites\"" >> ${OUTPATH_ROOT}/05_cand_sites/${UNIQ_ID}_05candSites.sh
echo -e "    /public/home/Songlab/xial/miniconda3/bin/python ~/scripts/pipeline/A2I_editing_pipe/A2I_pipe/vcf2RNAvcf.py -v \${OUTPATH_ROOT}/04_region/\${uniq_id}.ALU.filtered.bcf -d ${DEPTH} -r ALU -o \${OUTpath}/\${uniq_id}_Alu.all_candsites > \${OUTpath}/\${uniq_id}.ALU_sites.log 2>&1 &" >> ${OUTPATH_ROOT}/05_cand_sites/${UNIQ_ID}_05candSites.sh
echo -e "    /public/home/Songlab/xial/miniconda3/bin/python ~/scripts/pipeline/A2I_editing_pipe/A2I_pipe/vcf2RNAvcf.py -v \${OUTPATH_ROOT}/04_region/\${uniq_id}.nonALU_rep.filtered.bcf -d ${DEPTH} -r NONALU -o \${OUTpath}/\${uniq_id}_NoAlu_rep.all_candsites -R ${REF_FA} -sr ${simple_repeat} > \${OUTpath}/\${uniq_id}.nonALU_sites.log 2>&1 &" >> ${OUTPATH_ROOT}/05_cand_sites/${UNIQ_ID}_05candSites.sh
echo -e "    /public/home/Songlab/xial/miniconda3/bin/python ~/scripts/pipeline/A2I_editing_pipe/A2I_pipe/vcf2RNAvcf.py -v \${OUTPATH_ROOT}/04_region/\${uniq_id}.other.filtered.bcf -d ${DEPTH} -r NONALU -o \${OUTpath}/\${uniq_id}_other.all_candsites -R ${REF_FA} -sr ${simple_repeat} > \${OUTpath}/\${uniq_id}.other_sites.log 2>&1 &" >> ${OUTPATH_ROOT}/05_cand_sites/${UNIQ_ID}_05candSites.sh
echo -e "    wait" >> ${OUTPATH_ROOT}/05_cand_sites/${UNIQ_ID}_05candSites.sh
echo -e "}" >> ${OUTPATH_ROOT}/05_cand_sites/${UNIQ_ID}_05candSites.sh
echo -e "\n\n"
for FASTQ in ${INPATH}/*.fastq;
do
    uniq_id=${UNIQ_ID}_`basename ${FASTQ%.*}`
    OUTpath=${OUTPATH_ROOT}/05_cand_sites
    echo -e "candSites_process ${OUTPATH_ROOT} ${uniq_id} ${OUTpath} &" >> ${OUTPATH_ROOT}/05_cand_sites/${UNIQ_ID}_05candSites.sh
done
echo -e "wait" >> ${OUTPATH_ROOT}/05_cand_sites/${UNIQ_ID}_05candSites.sh
    
# sh ${OUTPATH_ROOT}/00_Nanoplot/${UNIQ_ID}_00Nanoplot.sh # > ${OUTPATH_ROOT}/00_Nanoplot/${UNIQ_ID}_00Nanoplot.log 2>&1
sh ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.sh # > ${OUTPATH_ROOT}/01_porechop/${UNIQ_ID}_01porechop.log 2>&1
sh ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.sh # > ${OUTPATH_ROOT}/02_Nanofilt/${UNIQ_ID}_02Nanofilt.log 2>&1
sh ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03minimap2.sh # > ${OUTPATH_ROOT}/03_minimap2/${UNIQ_ID}_03minimap2.log 2>&1
sh ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh # > ${OUTPATH_ROOT}/04_region/${UNIQ_ID}_04region.sh 2>&1
sh ${OUTPATH_ROOT}/05_cand_sites/${UNIQ_ID}_05candSites.sh
