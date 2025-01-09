library(ggplot2)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
inputFile <- args[1] # test_paramater.ALU.de_snpdb.candsite.txt
outPrefix <- as.numeric(args[2])
#prefix <- paste(strsplit(inputFile, "\\.")[[1]][1:length(strsplit(inputFile, "\\.")[[1]])-1], collapse='.')

gold.standard <- read.table('/public/home/Songlab/songyl/database/editing_sites/Merged/human/human_all_merged_editing-sites_add-strand', sep='\t', header=F)
colnames(gold.standard) <- c('CHROM', 'POS', 'STRAND')
gold.standard$ID <- paste(gold.standard$CHROM, gold.standard$POS, sep='_')

data.nano <- read.table(inputFile, header=T, sep='\t', row.names=NULL)
data.nano <- data.nano[data.nano$SEQ_ERROR_P > 0.05, ]
data.nano <- data.nano[data.nano$DELETION_NUM == 0 & data.nano$BP_RATIO > 5, ]
data.nano$MUT <- paste(toupper(data.nano$REF), toupper(data.nano$ALT), sep='>')
data.nano <- data.nano[(!(grepl('N', data.nano$MUT))),]
data.nano$REF_NUM <- as.numeric(sapply(data.nano$READS_NUM, function(x) strsplit(x, ",")[[1]][1]))
data.nano$ALT_NUM <- as.numeric(sapply(data.nano$READS_NUM, function(x) strsplit(x, ",")[[1]][2]))
#data.nano <- data.nano[data.nano$MUT %in% c('C>T', 'G>A'), ]
data.nano$ID <- paste(data.nano$CHROM, data.nano$POS, sep='_')
data.nano$GS <- 0
data.nano[data.nano$ID %in% gold.standard$ID, 'GS'] <- 1
data.nano$QUAL_pvalue <- sapply(data.nano$ALT_QUALITY, function(x) paste(sapply(as.numeric(strsplit(x, ",")[[1]]), function(x) round(10^((-1)*x/10),4)), collapse=','))
data.nano$QUAL_MEAN <- sapply(data.nano$QUAL_pvalue, function(x) mean(as.numeric(strsplit(x, ",")[[1]])))
true_pos <- data.nano[data.nano$GS == 1, 'ID']
false_pos <- data.nano[data.nano$GS == 0, 'ID']
print(paste(c('Total pos: ', nrow(data.nano)), collapse=''))
print(paste(c('TRUE pos: ', length(true_pos)), collapse=''))
print(paste(c('FALSE pos: ', length(false_pos)), collapse=''))
rand_num <- round(nrow(data.nano)/(10**(floor(log10(nrow(data.nano))))), digits=1) * 10**(floor(log10(nrow(data.nano)))-1)
print(paste(c('Random sample: ', rand_num), collapse=''))

rand_iter <- 100
QUAL_MEAN_iter <- c(); BP_RATIO_iter <- c(); HOMOPOLYMER_iter <- c(); DELETION_NUM_iter <- c(); PUZZLE_NUM_iter <- c(); SEQ_ERROR_P_iter <- c(); DELETION_NUM_RATIO_iter <- c(); PUZZLE_NUM_RATIO_iter <- c()
stored_ID_list <- c(); score_list <- c()
for (i in seq(1,rand_iter)){
    true_sample <- sample(1:length(true_pos), rand_num)
    false_sample <- sample(1:length(false_pos), rand_num)
    rand_true_pos <- true_pos[true_sample];
    rand_false_pos <- false_pos[false_sample];
    rand_test <- data.nano[data.nano$ID %in% c(rand_true_pos, rand_false_pos), c('EDITLEVEL', 'COVERAGE', 'COVERAGE_noSKIP', 'DELETION_NUM', 'BP_RATIO', 'SKIP', 'PUZZLE_NUM', 'SIMPLE_REPEAT', 'HOMOPOLYMER', 'SPLICE_NT', 'SEQ_ERROR_P', 'GS', 'QUAL_MEAN', 'ALT_NUM')]
    glm_test <- glm(GS~DELETION_NUM+BP_RATIO+PUZZLE_NUM+HOMOPOLYMER+SEQ_ERROR_P+QUAL_MEAN, data=rand_test, family=binomial(link='logit'))
    other_test_ID <- data.nano[(!(data.nano$ID %in% c(rand_true_pos, rand_false_pos))), 'ID']
    other_sample <- sample(1:length(other_test_ID), 10000)
    other_test <- data.nano[data.nano$ID %in% other_test_ID[other_sample], c('GS', 'DELETION_NUM', 'COVERAGE', 'BP_RATIO', 'PUZZLE_NUM', 'HOMOPOLYMER', 'SEQ_ERROR_P', 'ID', 'QUAL_MEAN')]
    other_test$score <- predict(glm_test, newdata=other_test, type='response')
    
    score_t <- fivenum(other_test$score)[4]
    filtered_other_test <- other_test[other_test$score > score_t, ]
    if (length(stored_ID_list) == 0){stored_ID_list <- filtered_other_test$ID
    } else {stored_ID_list <- c(stored_ID_list, filtered_other_test$ID)
    }
    if (length(score_list) == 0){score_list <- filtered_other_test$score
    } else {score_list <- c(score_list, filtered_other_test$score)}
    #output.data <- data.nano[data.nano$ID %in% stored_ID, ]
    #output.type <- aggregate(output.data$MUT, by=list(output.data$MUT), length)
    #colnames(output.type) <- c('MUT', 'number')
    #output.type$total <- length(stored_ID)
    #output.type$RATIO <- output.type$number/output.type$total

    QUAL_MEAN_iter <- c(QUAL_MEAN_iter, glm_test$coefficients['QUAL_MEAN'])
    BP_RATIO_iter <- c(BP_RATIO_iter, glm_test$coefficients['BP_RATIO'])
    HOMOPOLYMER_iter <- c(HOMOPOLYMER_iter, glm_test$coefficients['HOMOPOLYMER'])
    DELETION_NUM_iter <- c(DELETION_NUM_iter, glm_test$coefficients['DELETION_NUM'])
    PUZZLE_NUM_iter <- c(PUZZLE_NUM_iter, glm_test$coefficients['PUZZLE_NUM'])
    SEQ_ERROR_P_iter <- c(SEQ_ERROR_P_iter, glm_test$coefficients['SEQ_ERROR_P'])
    DELETION_NUM_RATIO_iter <- c(DELETION_NUM_RATIO_iter, glm_test$coefficients['DELETION_NUM:COVERAGE'])
    PUZZLE_NUM_RATIO_iter <- c(PUZZLE_NUM_RATIO_iter, glm_test$coefficients['COVERAGE:PUZZLE_NUM'])   
}

ITER.df <- data.frame(ITER=seq(1,rand_iter), QUAL_MEAN=QUAL_MEAN_iter, BP_RATIO=BP_RATIO_iter, HOMOPOLYMER=HOMOPOLYMER_iter, DELETION_NUM=DELETION_NUM_iter, PUZZLE_NUM=PUZZLE_NUM_iter, SEQ_ERROR_P=SEQ_ERROR_P_iter, DELETION_NUM_RATIO=DELETION_NUM_RATIO_iter, PUZZLE_NUM_RATIO=PUZZLE_NUM_RATIO_iter)
write.table(ITER.df, paste(c(prefix, 'paramater.txt'), collapse='.'), sep='\t', row.names=F, quote=F)
ITER.df <- read.table(paste(c(prefix, 'paramater.txt'), collapse='.'), sep='\t', header=T)
plot.df <- tidyr::pivot_longer(ITER.df, cols=2:ncol(ITER.df), values_to='Estimate', names_to='Paramater')
plot.df$X <- rep(1, nrow(plot.df))
p.para <- ggplot(plot.df, aes(x=X, y=Estimate)) + geom_boxplot() + egg::theme_article() + labs(x='', y='Estimate') +
    theme(axis.text = element_text(size=25), axis.title = element_text(size=30), strip.text=element_text(size=25), axis.ticks.length.y = unit(0.3, 'mm'), axis.text.x=element_blank()) +
    facet_wrap(.~Paramater, scales='free_y')
ggsave(paste(c(outPrefix, 'paramater.pdf'), collapse='.'), p.para, width=20, height=20)

pos_score.df <- data.frame(ID=stored_ID_list, SCORE=score_list)
pos_mean_score.df <- aggregate(pos_score.df$SCORE, by=list(pos_score.df$ID), mean)
colnames(pos_mean_score.df) <- c('ID', 'MEAN_SCORE')
pos_mean_score.df$MEAN_SCORE <- round(pos_mean_score.df$MEAN_SCORE, 4)
filtered.pos_score.df <- pos_mean_score.df[pos_mean_score.df$MEAN_SCORE > fivenum(pos_mean_score.df$MEAN_SCORE)[2], ]

output.type_func <- function(output.data){
    output.type <- aggregate(output.data$MUT, by=list(output.data$MUT), length)
    colnames(output.type) <- c('MUT', 'number')
    output.type$RATIO <- output.type$number / nrow(output.data)
    print(output.type)
}

output.data <- data.nano[data.nano$ID %in% filtered.pos_score.df$ID, ]
output.data.merge <- merge(output.data, filtered.pos_score.df, by='ID')
output.data.merge <- output.data.merge[, c(colnames(data.nano)[! colnames(data.nano) %in% c('ID', 'MUT', 'GS', 'QUAL_pvalue')], 'MEAN_SCORE')]
output.type <- aggregate(output.data$MUT, by=list(output.data$MUT), length)
colnames(output.type) <- c('MUT', 'number')
#output.type$total <- nrow(filtered.pos_score.df)
#output.type$RATIO <- output.type$number/output.type$total
output.type$RATIO <- output.type$number / nrow(output.data)
print(output.type)
write.table(output.data.merge, paste(c(outPrefix, 'filtered.txt'), collapse='.'), sep='\t', row.names=F, quote=F)

output.type[output.type$MUT %in% c('A>C', 'A>T', 'C>A', 'C>G', 'G>C', 'G>T', 'T>A', 'T>G'), 'MUT'] <- 'others'
output.type.new <- aggregate(output.type$number, by=list(output.type$MUT), sum)
colnames(output.type.new) <- c('MUT', 'SUM')
output.type.new$TOTAL <- sum(output.type.new$SUM)
output.type.new$RATIO <- output.type.new$SUM / output.type.new$TOTAL
output.type.new$MUT <- factor(output.type.new$MUT, levels=c('others', 'G>A', 'C>T', 'T>C', 'A>G'))
output.type.new$X <- rep('1', nrow(output.type.new))
p <- ggplot(output.type.new, aes(x=X, y=RATIO, fill=MUT)) + geom_bar(stat='identity', width=0.8) + scale_y_continuous(expand=expansion(mult = c(0,0))) +
    egg::theme_article() + scale_fill_manual(values=c('#D3D3D3', '#87CEFA', '#6495ED', '#FFDAB9', '#FF7F50')) + labs(x='', y='RATIO', title=paste(c('n=', unique(output.type.new$TOTAL)), collapse='')) +
    theme(axis.text.x=element_blank(), axis.text.y=element_text(size=25), axis.title.y=element_text(size=30), legend.title=element_blank(), legend.text=element_text(size=20), plot.title=element_text(size=25))
ggsave(paste(c(outPrefix, 'filtered.edit_type.pdf'), collapse='.'), p, width=6, height=10)

#mean_df <- c(mean(QUAL_MEAN_iter), mean(BP_RATIO_iter), mean(HOMOPOLYMER_iter), mean(DELETION_NUM_iter), mean(PUZZLE_NUM_iter), mean(SEQ_ERROR_P_iter), mean(DELETION_NUM_RATIO_iter), mean(PUZZLE_NUM_RATIO_iter))
#sd_df <- c(sd(QUAL_MEAN_iter), sd(BP_RATIO_iter), sd(HOMOPOLYMER_iter), sd(DELETION_NUM_iter), sd(PUZZLE_NUM_iter), sd(SEQ_ERROR_P_iter), sd(DELETION_NUM_RATIO_iter), sd(PUZZLE_NUM_RATIO_iter))
#para.df <- data.frame(PARAMATER=c('QUAL_MEAN', 'BP_RATIO', 'HOMOPOLYMER', 'DELETION_NUM', 'PUZZLE_NUM', 'SEQ_ERROR_P', 'DELETION_NUM:COVERAGE', 'PUZZLE_NUM:COVERAGE'), MEAN=mean_df, SD=sd_df)
#
#other_test_ID <- data.nano[(!(data.nano$ID %in% c(rand_true_pos, rand_false_pos))), 'ID']
#other_sample <- sample(1:length(other_test_ID), 10000)
#other_test <- data.nano[data.nano$ID %in% other_test_ID[other_sample], c('GS', 'DELETION_NUM', 'COVERAGE', 'BP_RATIO', 'PUZZLE_NUM', 'HOMOPOLYMER', 'SEQ_ERROR_P', 'ID', 'QUAL_MEAN')]
#other_test$score <- other_test$QUAL_MEAN * (-15.21)+other_test$BP_RATIO * (-6.813e-4)+other_test$HOMOPOLYMER*(-9.26e-1)+other_test$DELETION_NUM*(1.262e-1)+other_test$PUZZLE_NUM*(-3.726e-2)+other_test$SEQ_ERROR_P*2.322+(other_test$DELETION_NUM/other_test$COVERAGE)*(-1.01e-4) + other_test$PUZZLE_NUM/other_test$COVERAGE*(4.672e-5)
##other_test$score <- predict(glm_test, newdata=other_test, type='response')
#
#score_list <- c(); RATIO_list <- c()
#other_test$GS <- factor(other_test$GS, levels=c(0, 1))
##ggplot(other_test, aes(x=GS, y=score)) + geom_boxplot()
#for (score_t in seq(-0.3, 1.6, length.out=50)){
#    score_list <- c(score_list, score_t)
#    stored_ID <- other_test[other_test$score > score_t, 'ID']
#    output.data <- data.nano[data.nano$ID %in% stored_ID, ]
#    
#    output.type <- aggregate(output.data$MUT, by=list(output.data$MUT), length)
#    colnames(output.type) <- c('MUT', 'number')
#    output.type$total <- length(stored_ID)
#    output.type$RATIO <- output.type$number/output.type$total
#    RATIO_sum <- sum(output.type[output.type$MUT %in% c('A>G', 'T>C'), 'RATIO'])
#    RATIO_list <- c(RATIO_list, RATIO_sum)
#}
#ratio.df <- data.frame(SCORE=score_list, RATIO=RATIO_list)
#
#ggplot(ratio.df, aes(x=SCORE, y=RATIO)) + geom_point()
