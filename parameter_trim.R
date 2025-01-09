library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
inputFile <- args[1] # test_paramater.ALU.de_snpdb.candsite.txt
rand_sample_num <- as.numeric(args[2]) # randomly select several sites to generate the logistic model
prefix <- paste(strsplit(inputFile, "\\.")[[1]][1:length(strsplit(inputFile, "\\.")[[1]])-1], collapse='.')

gold.standard <- read.table('/public/home/Songlab/songyl/database/editing_sites/Merged/human/human_all_merged_editing-sites_add-strand', sep='\t', header=F)
colnames(gold.standard) <- c('CHROM', 'POS', 'STRAND')
gold.standard$ID <- paste(gold.standard$CHROM, gold.standard$POS, sep='_')

data.nano <- read.table(inputFile, header=T, sep='\t', row.names=NULL)
data.nano <- data.nano[data.nano$SEQ_ERROR_P > 0.05, ]
data.nano$MUT <- paste(toupper(data.nano$REF), toupper(data.nano$ALT), sep='>')
data.nano$ALT_NUM <- as.numeric(sapply(data.nano$READS_NUM, function(x) strsplit(x, ",")[[1]][2]))
#data.nano <- data.nano[data.nano$MUT %in% c('C>T', 'G>A'), ]
data.nano$ID <- paste(data.nano$CHROM, data.nano$POS, sep='_')
data.nano <- data.nano[data.nano$DELETION_NUM == 0 & data.nano$BP_RATIO > 5, ]
data.nano$GS <- 0
data.nano[data.nano$ID %in% gold.standard$ID, 'GS'] <- 1
data.nano$QUAL_pvalue <- sapply(data.nano$ALT_QUALITY, function(x) paste(sapply(as.numeric(strsplit(x, ",")[[1]]), function(x) round(10^((-1)*x/10),4)), collapse=','))
data.nano$QUAL_MEAN <- sapply(data.nano$QUAL_pvalue, function(x) mean(as.numeric(strsplit(x, ",")[[1]])))
data.nano$SKIP_NUM <- data.nano$COVERAGE - data.nano$COVERAGE_noSKIP
true_pos <- data.nano[data.nano$GS == 1, 'ID']
false_pos <- data.nano[data.nano$GS == 0, 'ID']

parameters <- c('SKIP_NUM', 'PUZZLE_NUM', 'HOMOPOLYMER', 'SEQ_ERROR_P', 'QUAL_MEAN', 'ALT_NUM', 'BP_RATIO')
QUAL_MEAN_iter <- c(); HOMOPOLYMER_iter <- c(); PUZZLE_NUM_iter <- c(); SEQ_ERROR_P_iter <- c(); ALT_NUM_iter <- c(); SKIP_NUM_iter <- c()
final_parameters <- c()
#for (i in seq(1,100)){
rand_true_pos <- sample(true_pos, rand_sample_num);
rand_false_pos <- sample(false_pos, rand_sample_num);
rand_test <- data.nano[data.nano$ID %in% c(rand_true_pos, rand_false_pos),]
glm_test <- glm(GS~DELETION_NUM+SKIP_NUM+BP_RATIO+PUZZLE_NUM+HOMOPOLYMER+SEQ_ERROR_P+QUAL_MEAN+ALT_NUM, data=rand_test, family=binomial(link='logit'))
final_parameters <- union(final_parameters, parameters[sapply(parameters, function(x) coef(summary(glm_test))[x, 4] < 0.05)])
    #parameters <- parameters[sapply(parameters, function(x) coef(summary(glm_test))[x, 4] < 0.05)]

#    QUAL_MEAN_iter <- c(QUAL_MEAN_iter, glm_test$coefficients['QUAL_MEAN'])
#    HOMOPOLYMER_iter <- c(HOMOPOLYMER_iter, glm_test$coefficients['HOMOPOLYMER'])
#    PUZZLE_NUM_iter <- c(PUZZLE_NUM_iter, glm_test$coefficients['PUZZLE_NUM'])
#    SEQ_ERROR_P_iter <- c(SEQ_ERROR_P_iter, glm_test$coefficients['SEQ_ERROR_P'])
#    ALT_NUM_iter <- c(ALT_NUM_iter, glm_test$coefficients['ALT_NUM'])
#    SKIP_NUM_iter <- c(SKIP_NUM_iter, glm_test$coefficients['SKIP_NUM'])
#}

#parameter.df <- data.frame(ITER=seq(1, length(QUAL_MEAN_iter)), QUAL_MEAN=QUAL_MEAN_iter, HOMOPOLYMER=HOMOPOLYMER_iter, PUZZLE_NUM=PUZZLE_NUM_iter, SEQ_ERROR_P=SEQ_ERROR_P_iter, ALT_NUM=ALT_NUM_iter, SKIP_NUM=SKIP_NUM_iter)
#final_parameter <- c()
#for (Feature in colnames(parameter.df)[2:ncol(paramater.df)]){
#    if (all(sapply(parameter.df[[Feature]], function(x) x<0)) || all(sapply(parameter.df[[Feature]], function(x) x>0))){
#        final_parameter <- c(final_parameter, Feature)
#    }
#}
#print(summary(glm_test))
#parameters <- c('DELETION_NUM', 'SKIP_NUM', 'BP_RATIO', 'PUZZLE_NUM', 'HOMOPOLYMER', 'SEQ_ERROR_P', 'QUAL_MEAN', 'ALT_NUM')
parameters <- c('SKIP_NUM', 'PUZZLE_NUM', 'HOMOPOLYMER', 'SEQ_ERROR_P', 'QUAL_MEAN', 'ALT_NUM')
parameters <- final_parameters
data.nano.filt <- data.nano
for (para in parameters) {
    #print(paste(c(para, mean(parameter.df[[para]])), collapse=','))
    #print(fivenum(rand_test[[para]]))
    #print(dim(data.nano.filt))
    if (coef(summary(glm_test))[para, 4] < 0.05) {
        if (glm_test$coefficients[[para]] > 0){
            #para_thre <- fivenum(rand_test[[para]])[2]
            para_thre <- fivenum(data.nano[[para]])[2]
            print(para_thre)
            data.nano.filt <- data.nano.filt[data.nano.filt[[para]] >= para_thre, ]
        } else {
            para_thre <- fivenum(data.nano[[para]])[4]
            print(para_thre)
            if (para_thre == 0){
                data.nano.filt <- data.nano.filt[data.nano.filt[[para]] == 0, ]
            } else {
                data.nano.filt <- data.nano.filt[data.nano.filt[[para]] <= para_thre, ]
            }
        }
    }
}

print(head(data.nano.filt))
output.type <- aggregate(data.nano.filt$MUT, by=list(data.nano.filt$MUT), length)
colnames(output.type) <- c('MUT', 'number')
output.type$total <- nrow(data.nano.filt)
output.type$RATIO <- output.type$number/output.type$total
print(output.type)

print(prefix)
write.table(data.nano.filt, paste(c(prefix, '.filtered.txt'), collapase=''), sep='\t', row.names=F)
