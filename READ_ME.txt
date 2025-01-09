Scripts for A2I_analysis
For analysing A2I RNA editing, including NGS A2I editing sites calling and level calling, and long-reads sequencing-based site calling, especially the Nanopore sequencing.

A2I_analysis.py:
    Python script for summarizing and filtering called editing sites
A2I_editing_pipeline.sh:
    Shell script for pipelining the calling processes according to Song et al processes
chisq_test.singleDatasets.R:
    R script for conducting chi-square test based on number of edited and unedited reads
downstream_summary_script:
    Shell-based pipeline for downstream plotting with called A2I sites
edit_level_replicate.R:
    R script for plotting the Venn-diagram, boxplot and ggpair plots to compare sites called within repeats
model_construction.py:
    Python script for filtering A2I sites using Long-reads sequencing technology with machine learning models
model_construction.R:
    R script for basically ploting distribution of paramaters generated in model_construction.py
pipe_levels.sh:
    Pipeline for NGS-based editing level calculation
