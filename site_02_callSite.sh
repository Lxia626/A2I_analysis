#! /usr/bin/env bash

#script used to pick up shared variants between two samples

cwd_path=`pwd`
filepath=$1 # absolute file path
Alu=$2 # {Alu, nonAlu, Both}
unique_id=$3
bam=${filepath}/${unique_id}.bam # merged file names with absolute path
sjmname=$bam".sjm"
if [[ ! -d ${filepath}/logs_site02 ]]; then mkdir -p ${filepath}/logs_site02; fi
if [[ -f ${bam} ]]; then rm ${bam}; fi
if [[ -f ${sjmname} ]]; then rm ${sjmname}; fi

#a=(`seq 1 22` 'X' 'Y' 'M')
a=(`seq 1 22` 'X' 'Y')
referencegenome="/public/home/Songlab/songyl/database/ref_genomes/Hg19/hg19_softmasked.fa"
hum_annotion="/public/home/Songlab/xial/scripts/pipeline_songyl/editing_call_pipeline/example_files/AllGeneAnno_hg19.txt"

echo -e "job_begin\n name ${unique_id}_MERGE\n memory 6G\n directory $filepath\n cmd_begin" >> ${sjmname}
bam_to_merge=""
for BAM in $filepath/starsplice_${unique_id}*/*_GATKRECAL.filtered.bam; do
    bam_to_merge=${bam_to_merge}" ${BAM}" >> ${sjmname}
done
echo -e "/public/usr/bin/samtools merge $bam ${bam_to_merge};" >> ${sjmname}
echo -e "/public/usr/bin/samtools index $bam;" >> ${sjmname}
echo -e "cmd_end\njob_end\n" >> ${sjmname}

SPLIT_ALU(){

    ALU=$1 # {Alu, nonAlu, Both}
    sjmname=$2
    candname=$3

    if [[ ${ALU} == "nonAlu" ]]; then
        suffix="NoAlu"
        echo -e "/public/usr/bin/bedtools intersect -a $candname.snp3.tmp -b /public/home/Songlab/songyl/database/hg19_Alu.bed -v > $candname.snp3.${suffix};" >> ${sjmname}
    else
        suffix="Alu"
        echo -e "/public/usr/bin/bedtools intersect -a $candname.snp3.tmp -b /public/home/Songlab/songyl/database/hg19_Alu.bed > $candname.snp3.${suffix};" >> ${sjmname}
    fi

    echo -e "sed s/\"\@\"/\"\\t\"/g $candname.snp3.${suffix} > $candname.snp3.${suffix}2;" >> ${sjmname}
    echo -e "awk -v OFS='\t' '{print \$1,\$3,\$4,\$5,\$6,\$7}' $candname.snp3.${suffix}2 > $candname.snp3.${suffix}2.format;" >> ${sjmname}
    echo -e "/usr/bin/perl /public/home/Songlab/xial/scripts/pipeline_songyl/editing_call_pipeline/Remove_mismatch_first6base.pl $candname.snp3.${suffix}2.format chr$key.filter.bam $candname.snp3.${suffix}2.format.r6;" >> ${sjmname}
    if [[ ${ALU} == "nonAlu" ]]; then
        echo -e "/usr/bin/perl /public/home/Songlab/xial/scripts/pipeline_songyl/editing_call_pipeline/filter_simple_repeats.pl $candname.snp3.${suffix}2.format.r6 /public/home/Songlab/songyl/database/hg19_SimpleRepeat.txt > $candname.snp3.${suffix}2.format.r6.rs;" >> ${sjmname}
        echo -e "/usr/bin/perl /public/home/Songlab/xial/scripts/pipeline_songyl/editing_call_pipeline/filter_splicing_artefacts_intronok_VR2.pl $candname.snp3.${suffix}2.format.r6.rs $hum_annotion > $candname.snp3.${suffix}2.format.r6.rs.rj;" >> ${sjmname}
        echo -e "/usr/bin/perl /public/home/Songlab/xial/scripts/pipeline_songyl/editing_call_pipeline/RemoveHomoNucleotides_forwei_5.pl $candname.snp3.${suffix}2.format.r6.rs.rj $candname.snp3.${suffix}2.format.r6.rs.rj.rho;" >> ${sjmname}
        echo -e "/usr/bin/perl /public/home/Songlab/xial/scripts/pipeline_songyl/editing_call_pipeline/Blat_alignedreads_onlychr_scorelimit_comparediscard_sanger.pl $candname.snp3.${suffix}2.format.r6.rs.rj.rho chr$key.filter.bam $candname.snp3.${suffix}2.format.r6.rs.rj.rho.blat;" >> ${sjmname}
        echo -e "/usr/bin/perl /public/home/Songlab/xial/scripts/pipeline_songyl/editing_call_pipeline/AnnotateGene_unstranded.pl $candname.snp3.${suffix}2.format.r6.rs.rj.rho.blat $hum_annotion > $candname.snp3.${suffix}2.format.r6.rs.rj.rho.blat_strand;" >> ${sjmname}
        echo -e "grep -v 'uc011aim.1' $candname.snp3.${suffix}2.format.r6.rs.rj.rho.blat_strand | grep -v 'uc010ytr.1' | grep -v 'uc010tyt.1' | grep -v 'uc010fhm.2' | grep -v 'uc002stl.2' | grep -v 'IG.V' > $candname.snp3.${suffix}2.format.r6.rs.rj.rho.blat_strand.rgene;" >> ${sjmname}
    else
        echo -e "/usr/bin/perl /public/home/Songlab/xial/scripts/pipeline_songyl/editing_call_pipeline/AnnotateGene_unstranded.pl $candname.snp3.${suffix}2.format.r6 $hum_annotion > $candname.snp3.${suffix}2.format.r6_strand;" >> ${sjmname}
        echo -e "grep -v 'uc011aim.1' $candname.snp3.${suffix}2.format.r6_strand | grep -v 'uc010ytr.1' | grep -v 'uc010tyt.1' | grep -v 'uc010fhm.2' | grep -v 'uc002stl.2' | grep -v 'IG.V' > $candname.snp3.${suffix}2.format.r6_strand.rgene;" >> ${sjmname}
    fi
    
}

for key in ${a[@]}; do

    pilename=$key".pileup"
    candname=$key".candsite.txt"

    echo -e "job_begin\n name ${unique_id}_SPLIT$key\n memory 30G\n  directory $filepath\n cmd_begin" >> ${sjmname}
    echo -e "/public/usr/bin/samtools view -bh $bam  chr$key -o chr$key.bam;"  >> ${sjmname}
    echo -e "/public/usr/bin/samtools view -b -F 1024 -q 20 -o chr$key.filter.bam chr$key.bam;" >> ${sjmname}
    echo -e "/public/usr/bin/samtools index chr$key.filter.bam;" >> ${sjmname}
    echo -e "/public/usr/bin/samtools mpileup -A -B -d 10000 -f $referencegenome -q 20 -Q 0 chr$key.filter.bam > $pilename;" >> ${sjmname}
    echo -e "/usr/bin/perl /public/home/Songlab/xial/scripts/pipeline_songyl/editing_call_pipeline/snp_caller2_sanger.pl $pilename > $candname;" >> ${sjmname} ###### CHANGE SNP CALLER PARAMETERS!!!!!!!!!!!!

    echo -e "/usr/bin/perl /public/home/Songlab/xial/scripts/pipeline_songyl/editing_call_pipeline/Compare_DARNED.pl $candname /public/home/Songlab/xial/Work/database/hg19/dbSNP151_SNV.site.txt $candname.snp1;" >> ${sjmname}
    echo -e "/usr/bin/perl /public/home/Songlab/xial/scripts/pipeline_songyl/editing_call_pipeline/Compare_DARNED.pl $candname.snp1 /public/home/Songlab/songyl/database/hg19_1000genomes_newSNPCalls.txt $candname.snp2;" >> ${sjmname}
    echo -e "/usr/bin/perl /public/home/Songlab/xial/scripts/pipeline_songyl/editing_call_pipeline/Compare_DARNED.pl $candname.snp2 /public/home/Songlab/songyl/database/UWash_variants.txt $candname.snp3;" >> ${sjmname}
    echo -e "awk -v OFS='\t' '{print \$1,\$2-1,\$2,\$3\"@\"\$4\"@\"\$5\"@\"\$6}' $candname.snp3 > $candname.snp3.tmp;" >> ${sjmname}

    if [[ ${Alu} == "Alu" ]]; then
        SPLIT_ALU Alu ${sjmname} ${candname}
    elif [[ ${Alu} == "nonAlu" ]]; then
        SPLIT_ALU nonAlu ${sjmname} ${candname}
    else # Both
        SPLIT_ALU Alu ${sjmname} ${candname}
        SPLIT_ALU nonAlu ${sjmname} ${candname}
    fi

    echo -e "cmd_end\njob_end\n" >> ${sjmname}

done

if [[ ${Alu} == "nonAlu" ]]; then
    echo -e "job_begin\n name ${unique_id}_CAT\n memory 3G\n directory $filepath\n cmd_begin" >> ${sjmname}
    echo -e "cat $filepath/*.blat_strand > ${unique_id}_all_candsites;" >> ${sjmname}
    echo -e "cmd_end\njob_end\n" >> ${sjmname}
elif [[ ${Alu} == "Alu" ]]; then
    echo -e "job_begin\n name ${unique_id}_CAT\n memory 3G\n directory $filepath\n cmd_begin" >> ${sjmname}
    echo -e "cat $filepath/*.r6_strand > ${unique_id}_all_candsites;" >> ${sjmname}
    echo -e "cmd_end\njob_end\n" >> ${sjmname}
else
    echo -e "job_begin\n name ${unique_id}_CAT_Alu\n memory 3G\n directory $filepath\n cmd_begin" >> ${sjmname}
    echo -e "cat $filepath/*.Alu*.r6_strand > ${unique_id}_Alu.all_candsites;" >> ${sjmname}
    echo -e "cmd_end\njob_end\n" >> ${sjmname}
    echo -e "job_begin\n name ${unique_id}_CAT_NoAlu\n memory 3G\n directory $filepath\n cmd_begin" >> ${sjmname}
    echo -e "cat $filepath/*.NoAlu*.blat_strand > ${unique_id}_NoAlu.all_candsites;" >> ${sjmname}
    echo -e "cmd_end\njob_end\n" >> ${sjmname}
fi

for key in ${a[@]}; do
    echo -e "order ${unique_id}_MERGE before ${unique_id}_SPLIT$key" >> ${sjmname}
done
if [[ ${Alu} == "Alu" || ${Alu} == "nonAlu" ]]; then
    for key in ${a[@]}; do
        echo -e "order ${unique_id}_SPLIT$key before ${unique_id}_CAT" >> ${sjmname}
    done
else
    for key in ${a[@]}; do
        echo -e "order ${unique_id}_SPLIT$key before ${unique_id}_CAT_Alu" >> ${sjmname}
    done
    echo -e "order ${unique_id}_CAT_Alu before ${unique_id}_CAT_NoAlu" >> ${sjmname}
fi
echo -e "log_dir $filepath/logs_site02\n" >> ${sjmname}
#sjm --max_running 10 ${sjmname}
#sjm ${sjmname}
