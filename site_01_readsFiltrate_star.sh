#! /usr/bin/env bash

filepath=$1 # filepath to current working directory
#user=$2 # user id when using qsub. {songyl, chenbj, xial}
stargenomeindex=$2 # INDEXed using STAR --runThreadN 10 --sjdbOverhang 99 --sjdbGTFfile gencode.v19.chr_patch_hapl_scaff.annotation.gtf --runMode genomeGenerate --genomeFastaFiles GRCh37.p13.genome.fa --genomeDir ./ according to Roth et al., Nat Methods, 2019

# divide each original FASTQ into several small part
for FASTQ in ${filepath}/*.fastq; do 
{
    fastqname=${FASTQ}
    newfile=${FASTQ}
    # basic file names
    fastqname=`basename ${FASTQ} | cut -d. -f1`
    foldername=${filepath}"/starsplice_"${fastqname}
    if [[ -d ${foldername} ]]; then rm -r ${foldername}; fi
    if [[ ! -d ${foldername} ]]; then mkdir ${foldername}; fi 
    if [[ ! -d ${foldername}/logs_site01 ]]; then 
        mkdir ${foldername}/logs_site01
    else
        rm ${foldername}/logs_site01/*
    fi
    if [[ ! -f $foldername"/splitfastq_aa" ]]; then 
        #qsub -cwd -b y -sync y split -l 50000000 $newfile $foldername/splitfastq # not to split FASTQ_
        ln -s ${FASTQ} $foldername"/splitfastq_aa"
    fi
}
done

# test whether qsub finished
wait
echo 'Starting Mapping...'

for FASTQ in ${filepath}/*.fastq; do

    fastqname=${FASTQ}
    newfile=${FASTQ}
    # basic file names
    fastqname=`basename ${FASTQ} | cut -d. -f1`
    foldername=${filepath}"/starsplice_"${fastqname}
    sortedname=${foldername}"/"${fastqname}"_sorted.bam"
    sjmname=${fastqname}"_star.sjm"
    logfile=${foldername}"/"${fastqname}"_star.log"
    if [[ -f ${sjmname} ]]; then rm ${sjmname}; fi
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

    # if [[ ! -f $foldername/splitfastq_aa ]]; then
    #     ln -s ${FASTQ} $foldername/splitfastq_aa
    # fi
    
    # initialize SJM script
    echo -e "job_begin\n name ${fastqname}_First\n memory 3G\n directory ${filepath}\n cmd_begin" >> ${sjmname}
    echo -e "/usr/bin/date >> $logfile;" >> ${sjmname}
    echo -e "cmd_end\njob_end\n" >> ${sjmname}

    for splitseqfile in $foldername/splitfastq_*; do

        tempname=${splitseqfile##*/}

        # mapping for each split FASTQs
        echo -e "job_begin\n name ${fastqname}_$tempname\n memory 40G\n directory $filepath\ncmd_begin" >> ${sjmname}
        #echo -e "/public/usr/bin/STAR --genomeDir $stargenomeindex --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 4 --readFilesIn $splitseqfile --outFileNamePrefix $splitseqfile"." --outTmpDir $foldername/$tempname"_tmpDir" > $splitseqfile.star.log;" >> ${sjmname}
        #echo -e "/public/usr/bin/STAR --genomeDir $stargenomeindex --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 4 --readFilesIn $splitseqfile --outFileNamePrefix $splitseqfile"." --outTmpDir $foldername/$tempname"_tmpDir" –alignIntronMax 1000000 –alignMatesGapMax 1000000 –alignSJoverhangMin 8 –outFilterMismatchNmax 999 –outFilterMismatchNoverReadLmax 0.1 –outFilterMultimapNmax 1 –outSAMattributes All > $splitseqfile.star.log;" >> ${sjmname} # according to Schaffer AA et al., NAR, 2020
        echo -e "/public/usr/bin/STAR --genomeDir $stargenomeindex --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 4 --readFilesIn $splitseqfile --outFileNamePrefix $splitseqfile"." --outTmpDir $foldername/$tempname"_tmpDir" --outFilterMatchNminOverLread 0.95 > $splitseqfile.star.log;" >> ${sjmname} # according to according to Roth et al., Nat Methods, 2019
        #echo -e "/public/usr/bin/STAR --genomeDir $stargenomeindex --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 4 --readFilesIn $splitseqfile --outFileNamePrefix $splitseqfile"." --outTmpDir $foldername/$tempname"_tmpDir" > $splitseqfile.star.log;" >> ${sjmname} # according to according to Mansi et al., NAR, 2020, followed by REDItools.
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

    ## correct the base quality and indel information
    #echo -e "/public/home/Songlab/xial/Work/software/gatk1/jdk-19.0.2/bin/java -jar /public/home/Songlab/xial/Work/software/gatk1/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $referencegenome -I $tmpname2 -o $intervalname;" >> ${sjmname}
    #echo -e "/public/home/Songlab/xial/Work/software/gatk1/jdk-19.0.2/bin/java -jar /public/home/Songlab/xial/Work/software/gatk1/GenomeAnalysisTK.jar -T IndelRealigner -R $referencegenome -I $tmpname2 -targetIntervals $intervalname -o $outfileRealigned;" >> ${sjmname}
    #echo -e "/public/home/Songlab/xial/Work/software/gatk1/jdk-19.0.2/bin/java -jar /public/home/Songlab/xial/Work/software/gatk1/GenomeAnalysisTK.jar -T CountCovariates -l INFO -D /public/home/Songlab/songyl/database/hg19_1000genomes_newSNPCalls.txt -I $outfileRealigned -R $referencegenome -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile $recalfile -dRG RG1 -dP illumina;" >> ${sjmname}
    #echo -e "/public/home/Songlab/xial/Work/software/gatk1/jdk-19.0.2/bin/java -jar /public/home/Songlab/xial/Work/software/gatk1/GenomeAnalysisTK.jar -T TableRecalibration -l INFO -I $outfileRealigned -R $referencegenome --out $outfileRecal -recalFile $recalfile -dRG RG1 -dP illumina;" >> ${sjmname}
    #echo -e "/public/usr/bin/samtools view -b -F 1024 -q 10 -o $bamtmp $outfileRecal;" >> ${sjmname}
    #echo -e "/usr/bin/date >> $logfile;" >> ${sjmname}
    #echo -e "cmd_end\njob_end\n" >> ${sjmname}

    # Next, using GATK3 for filtration
    #intervalname_3=${prefix}'_3.intervals'
    #outfileRealigned_3=${prefix}'_realignedGATK_3.bam'
    #outfileGrp_3=${prefix}'_recal_data_3.grp'
    #outfileRecal_3=${prefix}'_recalGATK_3.bam'
    #bamtmp_3=${prefix}'_GATKRECAL_3.filtered.bam'
    echo -e "/usr/bin/java -jar /public/home/Songlab/xial/Work/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $referencegenome -I $tmpname2 -o $intervalname -U ALLOW_N_CIGAR_READS -allowPotentiallyMisencodedQuals;" >> ${sjmname}
    echo -e "/usr/bin/java -jar /public/home/Songlab/xial/Work/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T IndelRealigner -R $referencegenome -I $tmpname2 -targetIntervals $intervalname -o $outfileRealigned -U ALLOW_N_CIGAR_READS -allowPotentiallyMisencodedQuals;" >> ${sjmname}
    echo -e "/usr/bin/java -jar /public/home/Songlab/xial/Work/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T BaseRecalibrator -R $referencegenome -I $outfileRealigned -knownSites /public/home/Songlab/songyl/database/GATK_ref/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf -knownSites /public/home/Songlab/songyl/database/GATK_ref/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -o $outfileGrp -U ALLOW_N_CIGAR_READS;" >> ${sjmname}
    echo -e "/usr/bin/java -jar /public/home/Songlab/xial/Work/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T PrintReads -R $referencegenome -I $outfileRealigned -BQSR $outfileGrp -o $outfileRecal -U ALLOW_N_CIGAR_READS;" >> ${sjmname}
    echo -e "/public/usr/bin/samtools view -b -F 1024 -q 10 -o $bamtmp $outfileRecal;" >> ${sjmname}
    echo -e "/usr/bin/date >> $logfile;" >> ${sjmname}
    echo -e "cmd_end\njob_end\n" >> ${sjmname}

    # total analysis order [What's the function?]
    for FILE in $foldername/splitfastq_*;do
        tempname=`basename ${FILE}`
        echo -e "order ${fastqname}_First before ${fastqname}_$tempname" >> ${sjmname}
        echo -e "order ${fastqname}_$tempname before ${fastqname}_FinishMapping" >> ${sjmname}
    done
    echo -e "order ${fastqname}_FinishMapping before ${fastqname}_MarkDuplicates" >> ${sjmname}
    echo -e "order ${fastqname}_MarkDuplicates before ${fastqname}_GATKCallEdits" >> ${sjmname}
    echo -e "log_dir $foldername/logs_site01" >> ${sjmname}
    # sjm $sjmname
done 
