**Scripts for A2I_analysis**  
For analysing A2I RNA editing, including NGS A2I editing sites calling and level calling, and long-reads sequencing-based site calling, especially the Nanopore sequencing.

`A2I_analysis.py`:  
	Python script for summarizing and filtering called editing sites  
`A2I_editing_pipeline.sh`:  
    Shell script for pipelining the calling processes according to Song et al processes  
`chisq_test.singleDatasets.R`:  
    R script for conducting chi-square test based on number of edited and unedited reads  
`downstream_summary_script`:  
    Shell-based pipeline for downstream plotting with called A2I sites  
`edit_level_replicate.R`:  
    R script for plotting the Venn-diagram, boxplot and ggpair plots to compare sites called within repeats  
`model_construction.py`:  
    Python script for filtering A2I sites using Long-reads sequencing technology with machine learning models  
`model_construction.R`:  
    R script for basically ploting distribution of paramaters generated in model_construction.py  
`mpileup_parser.py`:  
    Python script to parse the mpileup results generated by samtools mpileup  
`Nanopore_A2I_editing_pipe.sh`:  
    Shell-based pipeline for A2I sites analysis using Long-reads sequencing  
`Nanopore_editing_level.py`:  
    Python script for calculating editing level in Nanopore pipeline (maybe not used?)  
`parameter_trim*.R`:  
    Rscripts for filtering A2I sites called using Long-reads sequencing with different paramater trimming strategy  
`pipe_levels.sh`:  
    Pipeline for NGS-based editing level calculation  
`plot.R`:  
    Rscript for plotting mutant ratio in ALU, nonAlu_rep and others in single or multiple sample(s)  
`reditools.py`:  
    Python script to analyse REDItools results  
`site_01_readsFiltrate_star.sh and site_02_callSite.s`h:  
    Shell-based scripts to call sites using STAR, based on Song et al pipelines  
`STAR_FINAL_PIPELINE_all.sh, STAR_FINAL_PIPELINE.sh and sjm_combine.sh`:  
    Shell-based SJM file generator after using site_01_readsFiltrate_star.sh and site_02_callSite.sh. The two scripts share similar function  
`site_seq_err_correction.py`:  
    Python script to correct called sites accoring to reads number, based on Gabey et al., Nat Commun, 2022  
`vcf2RNAvcf.py`:  
    Python script to filter SNV sites generate by BCFtools (indexed)  
