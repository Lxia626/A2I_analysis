#! /usr/bin/env bash

filepath=$1 # filepath to current working directory
Alu=$2 # {Alu, nonAlu, Both}
unique_id=$3
stargenomeindex=$4 # INDEXed using STAR --runThreadN 10 --sjdbOverhang 99 --sjdbGTFfile gencode.v19.chr_patch_hapl_scaff.annotation.gtf --runMode genomeGenerate --genomeFastaFiles GRCh37.p13.genome.fa --genomeDir ./ according to Roth et al., Nat Methods, 2019
bam=${filepath}/${unique_id}.bam
sjmname=${filepath}/${unique_id}_total.sjm
#user=$2 # user id when using qsub. {songyl, xial}
    
# gtf: /public/home/Songlab/songyl/database/annotation/human/gencode.v19.chr_patch_hapl_scaff.annotation.gtf
# ref: /public/home/Songlab/songyl/database/genome/human/GencodeV19/STAR_index/with_GencodeV19_gtf/GRCh37.p13.genome.fa

# divide each original FASTQ into several small part
for FASTQ in ${filepath}/*.fastq; do 
{
    fastqname=${FASTQ}
    newfile=${FASTQ}
    # basic file names
    fastqname=`basename ${FASTQ} | cut -d. -f1`
    foldername=${filepath}"/starsplice_"${fastqname}
    if [[ ! -d ${foldername} ]]; then mkdir ${foldername}; fi 
    if [[ ! -d ${foldername}/logs_site01 ]]; then 
        mkdir ${foldername}/logs_site01
    else
        rm ${foldername}/logs_site01/*
    fi
    if [[ ! -f $foldername"/splitfastq_aa" ]]; then qsub -cwd -b y -sync y split -l 50000000 $newfile $foldername/splitfastq_; fi
}
done

# test whether qsub finished
wait

#######################
#### part1: mapping and reads filtration
#######################
if [[ -f ${sjmname} ]]; then rm ${sjmname}; fi
if [[ ! -d $filepath/logs_pipeline ]]; then 
    mkdir -p $filepath/logs_pipeline
else
    rm $filepath/logs_pipeline/*
fi

# initialize SJM script
echo -e "job_begin\n name ${unique_id}_First\n memory 3G\n directory ${filepath}\n cmd_begin" >> ${sjmname}
echo -e "/usr/bin/date;" >> ${sjmname}
echo -e "cmd_end\njob_end\n" >> ${sjmname}
a=1
for FASTQ in ${filepath}/*.fastq; do

    newfile=${FASTQ}
    # basic file names
    fastqname=`basename ${FASTQ} | cut -d. -f1`
    fastq_list[${a}]=${fastqname}
    a=`expr ${a} + 1`
    foldername=${filepath}"/starsplice_"${fastqname}
    sortedname=${foldername}"/"${fastqname}"_sorted.bam"
    logfile=${foldername}"/"${fastqname}"_star.log"
    if [[ ! -f ${logfile} ]]; then touch ${logfile}; fi

    # reference index
    # stargenomeindex="/public/home/Songlab/songyl/database/genome/human/GencodeV19/STAR_index/with_GencodeV19_gtf/"
    # stargenomeindex="/public/home/Songlab/xial/Work/Nanopore_RNA_editing/data/ylsong_GM12878/call_sites/Illumina/GM12878_new_call_sites/SRR3192400_lxia/call_site_star/STAR_index" # INDEXed using STAR --runThreadN 10 --sjdbOverhang 99 --sjdbGTFfile gencode.v19.chr_patch_hapl_scaff.annotation.gtf --runMode genomeGenerate --genomeFastaFiles GRCh37.p13.genome.fa --genomeDir ./ according to Roth et al., Nat Methods, 2019
    referencegenome="/public/home/Songlab/songyl/database/ref_genomes/Hg19/hg19_softmasked.fa"

    # output file names
    prefix=${foldername}"/"${fastqname}
    sortedprefix=${prefix}"_sorted"
    dedupname=${prefix}"_sorted.dedup.bam"
    tmpname=${prefix}"_sorted.bam.tmp"
    metricsname=${prefix}"_sorted.dedup.bam.dupmetrics"
    convertname=${prefix}"_converting.sam"
    convertname2=${prefix}"_converting2.sam"
    convertname3=${prefix}"_converting3.bam"
    dedupprefix=${prefix}"_sorted_final.dedup"
    bamtmp2=${prefix}"_sorted.dedup.filtered.bam"
    dedupnamefinal=${prefix}"_sorted_final.dedup.bam"
    intervalname=${prefix}".intervals"
    recalfile=${prefix}".recal_data.csv"
    outfileRealigned=${prefix}"_realignedGATK.bam"
    outfileRealignedIndex=${prefix}"_realignedGATK.bai"
    outfileRecal=${prefix}"_recalGATK.bam"
    outfileRecalIndex=${prefix}"_recalGATK.bai"
	tmpname2=${prefix}"_readgroupchange.bam"
    bamtmp=${prefix}"_GATKRECAL.filtered.bam"
    outfileGrp=${prefix}'_recal_data.grp'

    for splitseqfile in $foldername/splitfastq_*; do

        tempname=${splitseqfile##*/}

        # mapping for each split FASTQs
        echo -e "job_begin\n name ${fastqname}_${tempname}\n memory 40G\n directory $filepath\ncmd_begin" >> ${sjmname}
        #echo -e "/public/usr/bin/STAR --genomeDir $stargenomeindex --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 4 --readFilesIn $splitseqfile --outFileNamePrefix $splitseqfile"." --outTmpDir $foldername/$tempname"_tmpDir" > $splitseqfile.star.log;" >> ${sjmname}
        #echo -e "/public/usr/bin/STAR --genomeDir $stargenomeindex --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 4 --readFilesIn $splitseqfile --outFileNamePrefix $splitseqfile"." --outTmpDir $foldername/$tempname"_tmpDir" –alignIntronMax 1000000 –alignMatesGapMax 1000000 –alignSJoverhangMin 8 –outFilterMismatchNmax 999 –outFilterMismatchNoverReadLmax 0.1 –outFilterMultimapNmax 1 –outSAMattributes All > $splitseqfile.star.log;" >> ${sjmname} # according to Schaffer AA et al., NAR, 2020
        echo -e "/public/usr/bin/STAR --genomeDir $stargenomeindex --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 4 --readFilesIn $splitseqfile --outFileNamePrefix $splitseqfile"." --outTmpDir $foldername/$tempname"_tmpDir" --outFilterMatchNminOverLread 0.95 > $splitseqfile.star.log;" >> ${sjmname} # according to according to Roth et al., Nat Methods, 2019
        echo -e "cmd_end\njob_end\n" >> ${sjmname}

    done

    # merge splitted BAMs
    echo -e "job_begin\n name ${fastqname}_FinishMapping\n memory 3G\n directory $filepath\n cmd_begin" >> ${sjmname}
    if [[ -f $foldername/splitfastq_ab ]]; then
        to_merge=""
        for splitfastqx in $foldername/splitfastq_*; do
            to_merge="${to_merge} ${splitfastqx}.Aligned.sortedByCoord.out.bam"
        done
        echo -e "/public/usr/bin/samtools merge $sortedname $to_merge;" >> ${sjmname}
    else
        echo -e "mv $foldername/splitfastq_aa.Aligned.sortedByCoord.out.bam $sortedname;" >> ${sjmname}
    fi

    #echo -e "rm $foldername/splitfastq_*;" >> ${sjmname} # maybe commentted. In Step2, splitfastq_* are still useful.
    echo -e "/public/usr/bin/samtools flagstat $sortedname > $sortedname.flagstat;" >> ${sjmname}
    echo -e "cmd_end\njob_end\n" >> ${sjmname}

    # start analyze the merged BAMs
    echo -e "job_begin\n name ${fastqname}_MarkDuplicates\n memory 50G\n directory $filepath\n cmd_begin" >> ${sjmname}
    echo -e "/usr/bin/bash;\nexport _JAVA_OPTIONS=-Xmx2g;" >> ${sjmname}
    echo -e "/usr/bin/java -jar /public/home/Songlab/xial/bin/picard.jar AddOrReplaceReadGroups -I $sortedname -O $tmpname --SORT_ORDER coordinate --RGID RG1 --RGLB LB1 --RGPL illumina --RGSM SM1 --RGPU BC1 --VALIDATION_STRINGENCY LENIENT --TMP_DIR $foldername;" >> ${sjmname}
    echo -e "/usr/bin/java -jar /public/home/Songlab/xial/bin/picard.jar MarkDuplicates -I $tmpname -O $dedupname --METRICS_FILE $metricsname --VALIDATION_STRINGENCY LENIENT --TMP_DIR $foldername --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000;" >> ${sjmname}
    echo -e "cmd_end\njob_end\n" >> ${sjmname}

    # For editing calling, modify and filter the reads using picard and GATK
    # first of all, using GATK1 for filtration
    echo -e "job_begin\n name ${fastqname}_GATKCallEdits\n memory 50G\n directory $filepath\n cmd_begin" >> ${sjmname}
    echo -e "/public/usr/bin/samtools view -F 4 -o $convertname $dedupname;" >> ${sjmname}
    echo -e "/usr/bin/perl /public/home/Songlab/songyl/database/scripts/editing_call_pipeline/ConvertSpliceCoords.pl $convertname $convertname2;" >> ${sjmname}
    echo -e "/public/usr/bin/samtools view -bST $referencegenome -F 4 -o $convertname3 $convertname2;" >> ${sjmname}
    echo -e "/public/usr/bin/samtools sort -o $dedupnamefinal $convertname3;" >> ${sjmname}
    echo -e "/public/usr/bin/samtools view -hb -F 1024 -q 10 -o $bamtmp2 $dedupnamefinal;" >> ${sjmname}
    echo -e "/usr/bin/bash;\nexport _JAVA_OPTIONS=-Xmx2g;" >> ${sjmname}
    echo -e "/usr/bin/java -jar /public/home/Songlab/xial/bin/picard.jar AddOrReplaceReadGroups -I $bamtmp2 -O $tmpname2 --SORT_ORDER coordinate --RGID RG1 --RGLB LB1 --RGPL illumina --RGSM SM1 --RGPU BC1 --VALIDATION_STRINGENCY LENIENT --TMP_DIR $foldername;" >> ${sjmname}
    echo -e "/public/usr/bin/samtools index $tmpname2;" >> ${sjmname}

    echo -e "/usr/bin/java -jar /public/home/Songlab/xial/Work/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $referencegenome -I $tmpname2 -o $intervalname -U ALLOW_N_CIGAR_READS -allowPotentiallyMisencodedQuals;" >> ${sjmname}
    echo -e "/usr/bin/java -jar /public/home/Songlab/xial/Work/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T IndelRealigner -R $referencegenome -I $tmpname2 -targetIntervals $intervalname -o $outfileRealigned -U ALLOW_N_CIGAR_READS -allowPotentiallyMisencodedQuals;" >> ${sjmname}
    echo -e "/usr/bin/java -jar /public/home/Songlab/xial/Work/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T BaseRecalibrator -R $referencegenome -I $outfileRealigned -knownSites /public/home/Songlab/songyl/database/GATK_ref/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf -knownSites /public/home/Songlab/songyl/database/GATK_ref/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -o $outfileGrp -U ALLOW_N_CIGAR_READS;" >> ${sjmname}
    echo -e "/usr/bin/java -jar /public/home/Songlab/xial/Work/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T PrintReads -R $referencegenome -I $outfileRealigned -BQSR $outfileGrp -o $outfileRecal -U ALLOW_N_CIGAR_READS;" >> ${sjmname}
    echo -e "/public/usr/bin/samtools view -b -F 1024 -q 10 -o $bamtmp $outfileRecal;" >> ${sjmname}
    echo -e "/usr/bin/date >> $logfile;" >> ${sjmname}
    echo -e "cmd_end\njob_end\n" >> ${sjmname}

done 

#########################
##### part2: candidate sites filtration
#########################

echo -e "job_begin\n name ${unique_id}_copy_BAM\n memory 3G\n directory $filepath\n cmd_begin" >> ${sjmname}
for FASTQ in ${filepath}/*.fastq; do
    fastqname=`basename ${FASTQ} | cut -d. -f1`
    prefix=${filepath}"/starsplice_"${fastqname}"/"${fastqname}
    echo -e "ln -s ${prefix}_GATKRECAL.filtered.bam ${filepath};" >> ${sjmname}
done
echo -e "cmd_end\njob_end\n" >> ${sjmname}

if [[ -f ${bam} ]]; then rm ${bam}; fi

chrom_list=(`seq 1 22` 'X' 'Y' 'M')
referencegenome="/public/home/Songlab/songyl/database/ref_genomes/Hg19/hg19_softmasked.fa"
hum_annotion="/public/home/Songlab/xial/scripts/pipeline_songyl/editing_call_pipeline/example_files/AllGeneAnno_hg19.txt"

echo -e "job_begin\n name ${unique_id}_MERGE\n memory 6G\n directory $filepath\n cmd_begin" >> ${sjmname}
bam_to_merge=""
for FASTQ in ${filepath}/*.fastq; do    
    fastqname=`basename ${FASTQ} | cut -d. -f1`
    prefix=${filepath}"/"${fastqname}
    bam_to_merge=${bam_to_merge}" ${prefix}_GATKRECAL.filtered.bam" >> ${sjmname}
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
    echo -e "/usr/bin/perl /public/home/Songlab/xial/scripts/pipeline_songyl/editing_call_pipeline/filter_simple_repeats.pl $candname.snp3.${suffix}2.format.r6 /public/home/Songlab/songyl/database/hg19_SimpleRepeat.txt > $candname.snp3.${suffix}2.format.r6.rs;" >> ${sjmname}
    echo -e "/usr/bin/perl /public/home/Songlab/xial/scripts/pipeline_songyl/editing_call_pipeline/filter_splicing_artefacts_intronok_VR2.pl $candname.snp3.${suffix}2.format.r6.rs $hum_annotion > $candname.snp3.${suffix}2.format.r6.rs.rj;" >> ${sjmname}
    echo -e "/usr/bin/perl /public/home/Songlab/xial/scripts/pipeline_songyl/editing_call_pipeline/RemoveHomoNucleotides_forwei_5.pl $candname.snp3.${suffix}2.format.r6.rs.rj $candname.snp3.${suffix}2.format.r6.rs.rj.rho;" >> ${sjmname}
    echo -e "/usr/bin/perl /public/home/Songlab/xial/scripts/pipeline_songyl/editing_call_pipeline/Blat_alignedreads_onlychr_scorelimit_comparediscard_sanger.pl $candname.snp3.${suffix}2.format.r6.rs.rj.rho chr$key.filter.bam $candname.snp3.${suffix}2.format.r6.rs.rj.rho.blat;" >> ${sjmname}
    echo -e "/usr/bin/perl /public/home/Songlab/xial/scripts/pipeline_songyl/editing_call_pipeline/AnnotateGene_unstranded.pl $candname.snp3.${suffix}2.format.r6.rs.rj.rho.blat $hum_annotion > $candname.snp3.${suffix}2.format.r6.rs.rj.rho.blat_strand;" >> ${sjmname}
    echo -e "grep -v 'uc011aim.1' $candname.snp3.${suffix}2.format.r6.rs.rj.rho.blat_strand | grep -v 'uc010ytr.1' | grep -v 'uc010tyt.1' | grep -v 'uc010fhm.2' | grep -v 'uc002stl.2' | grep -v 'IG.V' > $candname.snp3.${suffix}2.format.r6.rs.rj.rho.blat_strand.rgene;" >> ${sjmname}

}

for key in ${chrom_list[@]}; do

    pilename=$key".pileup"
    candname=$key".candsite.txt"

    echo -e "job_begin\n name ${unique_id}_SPLIT$key\n memory 30G\n  directory $filepath\n cmd_begin" >> ${sjmname}
    echo -e "/public/usr/bin/samtools view -bh $bam  chr$key -o chr$key.bam;"  >> ${sjmname}
    echo -e "/public/usr/bin/samtools view -b -F 1024 -q 20 -o chr$key.filter.bam chr$key.bam;" >> ${sjmname}
    echo -e "/public/usr/bin/samtools index chr$key.filter.bam;" >> ${sjmname}
    echo -e "/public/usr/bin/samtools mpileup -A -B -d 10000 -f $referencegenome -q 20 -Q 0 chr$key.filter.bam > $pilename;" >> ${sjmname}
    echo -e "/usr/bin/perl /public/home/Songlab/xial/scripts/pipeline_songyl/editing_call_pipeline/snp_caller2_sanger.pl $pilename > $candname;" >> ${sjmname} ###### CHANGE SNP CALLER PARAMETERS!!!!!!!!!!!!

    echo -e "/usr/bin/perl /public/home/Songlab/xial/scripts/pipeline_songyl/editing_call_pipeline/Compare_DARNED.pl $candname /public/home/Songlab/xial/Work/database/dbSNP151_SNV.site.txt $candname.snp1;" >> ${sjmname}
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

if [[ ${Alu} == "Alu" || ${Alu} == "nonAlu" ]]; then
    echo -e "job_begin\n name ${unique_id}_CAT\n memory 3G\n directory $filepath\n cmd_begin" >> ${sjmname}
    echo -e "cat $filepath/*.blat_strand > all_candsites;" >> ${sjmname}
    echo -e "cmd_end\njob_end\n" >> ${sjmname}
else
    echo -e "job_begin\n name ${unique_id}_CAT_Alu\n memory 3G\n directory $filepath\n cmd_begin" >> ${sjmname}
    echo -e "cat $filepath/*.Alu*.blat_strand > Alu.all_candsites;" >> ${sjmname}
    echo -e "cmd_end\njob_end\n" >> ${sjmname}
    echo -e "job_begin\n name ${unique_id}_CAT_NoAlu\n memory 3G\n directory $filepath\n cmd_begin" >> ${sjmname}
    echo -e "cat $filepath/*.NoAlu*.blat_strand > NoAlu.all_candsites;" >> ${sjmname}
    echo -e "cmd_end\njob_end\n" >> ${sjmname}
fi

# whole pipeline
### part 1: mapping and filtration
for FASTQ_name in ${fastq_list[@]}; do
    for FILE in $foldername/splitfastq_*; do
        tempname=`basename ${FILE}`
        echo -e "order ${unique_id}_First before ${FASTQ_name}_$tempname" >> ${sjmname}
        echo -e "order ${FASTQ_name}_$tempname before ${FASTQ_name}_FinishMapping" >> ${sjmname}
    done
    echo -e "order ${FASTQ_name}_FinishMapping before ${FASTQ_name}_MarkDuplicates" >> ${sjmname}
    echo -e "order ${FASTQ_name}_MarkDuplicates before ${FASTQ_name}_GATKCallEdits" >> ${sjmname}
    echo -e "order ${FASTQ_name}_GATKCallEdits before ${unique_id}_copy_BAM" >> ${sjmname}
done

### part 2: site filtration
echo -e "order ${unique_id}_copy_BAM before ${unique_id}_MERGE" >> ${sjmname}
for key in ${chrom_list[@]}; do
    echo -e "order ${unique_id}_MERGE before ${unique_id}_SPLIT$key" >> ${sjmname}
done
if [[ ${Alu} == "Alu" || ${Alu} == "nonAlu" ]]; then
    for key in ${chrom_list[@]}; do
        echo -e "order ${unique_id}_SPLIT$key before ${unique_id}_CAT" >> ${sjmname}
    done
else
    for key in ${chrom_list[@]}; do
        echo -e "order ${unique_id}_SPLIT$key before ${unique_id}_CAT_Alu" >> ${sjmname}
    done
    echo -e "order ${unique_id}_CAT_Alu before ${unique_id}_CAT_NoAlu" >> ${sjmname}
fi
echo -e "log_dir $filepath/logs_pipeline\n" >> ${sjmname}
# sjm --max_running 10 ${sjmname}