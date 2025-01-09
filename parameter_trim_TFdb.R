library(ggplot2)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
inputFile <- args[1] # test_paramater.ALU.de_snpdb.candsite.txt
outPrefix <- as.numeric(args[2])
#prefix <- paste(strsplit(inputFile, "\\.")[[1]][1:length(strsplit(inputFile, "\\.")[[1]])-1], collapse='.')

#rna.standard <- read.table('/public/home/Songlab/songyl/database/editing_sites/Merged/human/human_all_merged_editing-sites_add-strand', sep='\t', header=F)
#colnames(rna.standard) <- c('CHROM', 'POS', 'STRAND')
#rna.standard$ID <- paste(rna.standard$CHROM, rna.standard$POS, sep='_')
#dna.standard <- read.table('/public/home/Songlab/xial/Work/database/hg19/SNP_all_dbSNP_1000Genomes_UWash.txt', sep='\t', header=F)
#colnames(dna.standard) <- c('CHROM', 'POS')
#dna.standard$ID <- paste(dna.standard$CHROM, dna.standard$POS, sep='_')

# candidate site annotation
data.nano <- read.table(inputFile, header=T, sep='\t', row.names=NULL)
data.nano <- data.nano[data.nano$SEQ_ERROR_P > 0.05, ]
data.nano <- data.nano[data.nano$DELETION_NUM == 0 & data.nano$BP_RATIO > 5 & data.nano$HOMOPOLYMER == 0 & data.nano$SIMPLE_REPEAT == 0, ]
data.nano$MUT <- paste(toupper(data.nano$REF), toupper(data.nano$ALT), sep='>')
data.nano <- data.nano[(!(grepl('N', data.nano$MUT))),]
data.nano$REF_NUM <- as.numeric(sapply(data.nano$READS_NUM, function(x) strsplit(x, ",")[[1]][1]))
data.nano$ALT_NUM <- as.numeric(sapply(data.nano$READS_NUM, function(x) strsplit(x, ",")[[1]][2]))
#data.nano <- data.nano[data.nano$MUT %in% c('C>T', 'G>A'), ]
data.nano$ID <- paste(data.nano$CHROM, data.nano$POS, sep='_')
data.nano$QUAL_pvalue <- sapply(data.nano$ALT_QUALITY, function(x) paste(sapply(as.numeric(strsplit(x, ",")[[1]]), function(x) round(10^((-1)*x/10),4)), collapse=','))
data.nano$SKIP_NUM <- data.nano$COVERAGE - data.nano$COVERAGE_noSKIP
data.nano$QUAL_MEAN <- sapply(data.nano$QUAL_pvalue, function(x) mean(as.numeric(strsplit(x, ",")[[1]])))
true_pos <- na.omit(data.nano[data.nano$DataBase == 'dbRNA', 'ID'])
false_pos <- na.omit(data.nano[data.nano$DataBase == 'dbSNP', 'ID'])
print(paste(c('Total pos: ', nrow(data.nano)), collapse=''))
print(paste(c('TRUE pos: ', length(true_pos)), collapse=''))
print(paste(c('FALSE pos: ', length(false_pos)), collapse=''))
rand_num <- round(nrow(data.nano)/(10**(floor(log10(nrow(data.nano))))), digits=1) * 10**(floor(log10(nrow(data.nano)))-1)
print(paste(c('Random sample: ', rand_num), collapse=''))

# model building and testing
training_set <- data.nano[data.nano$ID %in% c(true_pos, false_pos), ]
training_set$GS <- 0
training_set[training_set$ID %in% true_pos, 'GS'] <- 1
#glm_test <- glm(GS~PUZZLE_NUM+HOMOPOLYMER+SEQ_ERROR_P+QUAL_MEAN+SKIP_NUM+ALT_NUM, data=training_set, family=binomial(link='logit'))
glm_test <- glm(GS~PUZZLE_NUM/COVERAGE_noSKIP+HOMOPOLYMER+SEQ_ERROR_P+QUAL_MEAN+SKIP_NUM+ALT_NUM, data=training_set, family=binomial(link='logit'))
predit_df <- data.nano[(!(data.nano$ID %in% c(true_pos, false_pos))), ]
predit_df$score <- predict(glm_test, newdata=predit_df, type='response')

# site filtration
score_t <- fivenum(predit_df$score)[4]
filtered_ID <- predit_df[predit_df$score > score_t, 'ID']
filtered_df <- data.nano[data.nano$ID %in% filtered_ID, ]
out_predict_df <- rbind(filtered_df, data.nano[data.nano$ID %in% true_pos, ])
write.table(out_predict_df, paste(c(outPrefix, 'filtered.txt'), collapse='.'), sep='\t', row.names=F, quote=F)

# ratio of MUT
output.type <- aggregate(out_predict_df$MUT, by=list(out_predict_df$MUT), length)
colnames(output.type) <- c('MUT', 'number')
output.type$RATIO <- output.type$number / nrow(out_predict_df)
print(output.type)
dbRNA.type <- aggregate(na.omit(out_predict_df[out_predict_df$DataBase == 'dbRNA', 'MUT']), by=list(na.omit(out_predict_df[out_predict_df$DataBase == 'dbRNA', 'MUT'])), length)
colnames(dbRNA.type) <- c('MUT', 'number')
dbRNA.type$RATIO <- dbRNA.type$number / length(na.omit(out_predict_df[out_predict_df$DataBase == 'dbRNA', 'MUT']))
print(dbRNA.type)
predict.type <- aggregate(out_predict_df[is.na(out_predict_df$DataBase), 'MUT'], by=list(out_predict_df[is.na(out_predict_df$DataBase), 'MUT']), length)
colnames(predict.type) <- c('MUT', 'number')
predict.type$RATIO <- predict.type$number / length(out_predict_df[is.na(out_predict_df$DataBase), 'MUT'])
print(predict.type)

write.table(output.type, paste(c(outPrefix, 'filtered.edt_type.txt'), collapse='.'), sep='\t', row.names=F, quote=F)

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
