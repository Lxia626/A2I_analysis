#! /usr/bin/env Rscript

library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
inputFile <- args[1] # test_paramater.ALU.de_snpdb.candsite.txt
prefix <- args[2]
dir.create('./plot')
iter_num <- 100
data_feature <- read.table(inputFile, sep="\t", header=T, quote='')
data_feature$REF <- toupper(data_feature$REF)
data_feature$ID <- paste(data_feature$CHROM, data_feature$POS, sep='_')
data_feature$MUT <- paste(data_feature$REF, data_feature$ALT, sep='_')

# MOTIF based paramater summarizing
motif_df <- data_feature[data_feature$MUT %in% c('A_G','T_C'), 
    c('MOTIF', 'QUAL_MEAN_1','QUAL_STD_1','DEL_NUM_1','DEL_RATE_1','DEL_LEN_1','INS_NUM_1','INS_RATE_1','INS_LEN_1','DEL_SITE_NUM_1','DEL_SITE_RATIO_1','DEPTH_1','BP_1','SKIP_1','PUZZLE_NUM_1','SR_1','HOMO_1','SPLICE_1',
    'QUAL_MEAN_mid','QUAL_STD_mid','DEL_NUM_mid','DEL_RATE_mid','DEL_LEN_mid','INS_NUM_mid','INS_RATE_mid','INS_LEN_mid','DEL_SITE_NUM_mid','DEL_SITE_RATIO_mid','DEPTH_mid','BP_mid','SKIP_mid','PUZZLE_NUM_mid','SR_mid','HOMO_mid','SPLICE_mid',
    'QUAL_MEAN_3','QUAL_STD_3','DEL_NUM_3','DEL_RATE_3','DEL_LEN_3','INS_NUM_3','INS_RATE_3','INS_LEN_3','DEL_SITE_NUM_3','DEL_SITE_RATIO_3','DEPTH_3','BP_3','SKIP_3','PUZZLE_NUM_3','SR_3','HOMO_3','SPLICE_3')]
motif_df[, 2:ncol(motif_df)] <- apply(motif_df[, 2:ncol(motif_df)], 2, as.numeric)
motif_df$CLASS <- 0
motif_df[motif_df$MOTIF %in% c('TAG', 'CTA'), 'CLASS' ] <- 1
motif_df$MOTIF <- factor(motif_df$MOTIF, levels=c('TAG','CTA','TAC','GTA','TAA','TTA','TAT','ATA',
    'CAG','CTG','CAT','ATG','CAC','GTG','CAA','TTG',
    'GAG','CTC','GAC','GTC','GAA','TTC','GAT','ATC',
    'AAG','CTT','AAC','GTT','AAT','ATT','AAA','TTT'))
motif_df.plot <- tidyr::pivot_longer(motif_df, cols=2:ncol(motif_df), names_to='Paramater', values_to='Values')
p.paramater <- ggplot(motif_df.plot, aes(x=MOTIF, y=Values)) + geom_boxplot() + egg::theme_article() +
    theme(axis.text=element_text(size=20), axis.title=element_text(size=25), axis.ticks.length=unit(3,'mm'), axis.text.x=element_text(angle=90)) +
    facet_wrap(.~Paramater, scales='free_y') +
    theme(strip.text=element_text(size=15))
ggsave('./plot/Paramater_distribution.pdf', p.paramater, width=60, height=60, limitsize = FALSE)
p.paramater.no <- ggplot(motif_df.plot, aes(x=MOTIF, y=Values)) + geom_boxplot(outlier.shape=NA) + egg::theme_article() + 
    theme(axis.text=element_text(size=20), axis.title=element_text(size=25), axis.ticks.length=unit(3,'mm'), axis.text.x=element_text(angle=90)) +
    facet_grid(Paramater~., scales='free_y', space='free_y', shrink=TRUE) +
    theme(strip.text=element_text(size=15))
ggsave('./plot/Paramater_distribution.noOut.pdf', p.paramater.no, width=60, height=60, limitsize = FALSE)

# UAG as TRUE, others as FALSE
rand_num <- round(nrow(motif_df)*0.1, 0)
for (i in seq(1:iter_num)){
    true_sample <- sample(1:dim(motif_df[motif_df$CLASS == 1,])[1], rand_num)
    false_sample <- sample(1:dim(motif_df[motif_df$CLASS == 0,])[1], rand_num)
    rand_test <- rbind(true_sample, false_sample)
    glm_test <- glm(CLASS~QUAL_MEAN_1+QUAL_STD_1+DEL_NUM_1+DEL_RATE_1+DEL_LEN_1+INS_NUM_1+INS_RATE_1+INS_LEN_1+DEL_SITE_NUM_1+DEL_SITE_RATIO_1+DEPTH_1+BP_1+SKIP_1+PUZZLE_NUM_1+SR_1+HOMO_1+SPLICE_1+
        QUAL_MEAN_mid+QUAL_STD_mid+DEL_NUM_mid+DEL_RATE_mid+DEL_LEN_mid+INS_NUM_mid+INS_RATE_mid+INS_LEN_mid+DEL_SITE_NUM_mid+DEL_SITE_RATIO_mid+DEPTH_mid+BP_mid+SKIP_mid+PUZZLE_NUM_mid+SR_mid+HOMO_mid+SPLICE_mid+
        QUAL_MEAN_3+QUAL_STD_3+DEL_NUM_3+DEL_RATE_3+DEL_LEN_3+INS_NUM_3+INS_RATE_3+INS_LEN_3+DEL_SITE_NUM_3+DEL_SITE_RATIO_3+DEPTH_3+BP_3+SKIP_3+PUZZLE_NUM_3+SR_3+HOMO_3+SPLICE_3, data=rand_test, family=binomial(link='logit'))
}

# known sites as TRUE, other sites as FALSE
known_sites <- read.table('/public/home/Songlab/songyl/database/editing_sites/Merged/human/human_all_merged_editing-sites_add-strand', sep='\t', header=F)
colnames(known_sites) <- c('CHROM', 'POS', 'STRAND')
known_sites$ID <- paste(known_sites$CHROM, known_sites$POS, sep='_')
data_feature$KNOWN <- 0
data_feature[data_feature$ID %in% known_sites$ID, 'KNOWN'] <- 1
true_df <- data_feature[data_feature$KNOWN == 1,]
false_df <- data_feature[data_feature$KNOWN == 0,]

rand_num <- round(nrow(data_feature)*0.1, 0)
for (i in seq(1:iter_num)){
    true_sample <- sample(1:dim(true_df)[1], rand_num)
    false_sample <- sample(1:dim(false_df)[1], rand_num)
    rand_true <- true_df[true_sample, ]
    rand_false <- false_df[true_sample, ]
    rand_test <- rbind(rand_false, rand_true)
    glm_test <- glm(KNOWN~QUAL_MEAN_1+QUAL_STD_1+DEL_NUM_1+DEL_RATE_1+DEL_LEN_1+INS_NUM_1+INS_RATE_1+INS_LEN_1+DEL_SITE_NUM_1+DEL_SITE_RATIO_1+DEPTH_1+BP_1+SKIP_1+PUZZLE_NUM_1+SR_1+HOMO_1+SPLICE_1+
        QUAL_MEAN_mid+QUAL_STD_mid+DEL_NUM_mid+DEL_RATE_mid+DEL_LEN_mid+INS_NUM_mid+INS_RATE_mid+INS_LEN_mid+DEL_SITE_NUM_mid+DEL_SITE_RATIO_mid+DEPTH_mid+BP_mid+SKIP_mid+PUZZLE_NUM_mid+SR_mid+HOMO_mid+SPLICE_mid+
        QUAL_MEAN_3+QUAL_STD_3+DEL_NUM_3+DEL_RATE_3+DEL_LEN_3+INS_NUM_3+INS_RATE_3+INS_LEN_3+DEL_SITE_NUM_3+DEL_SITE_RATIO_3+DEPTH_3+BP_3+SKIP_3+PUZZLE_NUM_3+SR_3+HOMO_3+SPLICE_3, data=rand_test, family=binomial(link='logit'))
    

