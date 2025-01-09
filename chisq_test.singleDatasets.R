library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
level_dir <- args[1]
out_prefix <- args[2]

level_df.out <- read.table(paste(c('replicates_plot', paste(c(out_prefix,'filtered.edit_level.txt'),collapse='.')), collapse='/'), header=T, sep='\t')
sample_list <- sapply(list.files('tmp_dir/'), function(x) strsplit(x,'\\.')[[1]][1])

# sample level
data_region <- data.frame()
for (FILE in list.files('tmp_dir/')){
    tmp.df <- read.table(paste(c('tmp_dir/',FILE), collapse=''), header=F, sep='\t')
    tmp.df$V4 <- gsub('N/A', NA, tmp.df$V4)
    tmp.df <- tmp.df[tmp.df$V2 > 30, ]
    tmp.df$V4 <- as.numeric(tmp.df$V4)
    tmp.df <- tmp.df[tmp.df$V1 %in% level_df.out$ID, ]
    colnames(tmp.df) <- c('ID', 'COVERAGE', 'EDITREADS', 'EDITLEVEL', 'REGION')
    tmp.df$SAMPLE <- strsplit(FILE, '\\.')[[1]][1]
    data_region <- rbind(data_region, tmp.df[, c('ID', 'COVERAGE', 'EDITREADS', 'EDITLEVEL', 'SAMPLE')])
}
write.table(data_region, paste(c('replicates_plot', paste(c(out_prefix, 'data_region.txt'), collapse='.')), collapse='/'), quote=F, sep='\t', row.names=F)

comparison_list <- combn(sample_list, 2, simplify=F)
#data_region <- read.table(paste(c('replicates_plot', paste(c(out_prefix, 'data_region.txt'), collapse='.')), collapse='/'), header=T, sep='\t')
#ggplot(data_region, aes(x=SAMPLE, y=EDITLEVEL)) + geom_boxplot() + egg::theme_article() +
#    scale_y_continuous(limits=c(0, 0.3))
    
### dataset level
#treatment.df <- data.frame()
#for (TREATMENT in treatment){
#    treat.df <- data.frame()
#    for (i in seq(1,3)){
#        tmp.df <- read.table(paste(c('tmp_dir/', paste(c(paste(c(TREATMENT,i),collapse='_'), 'edit_level.txt'), collapse='.')), collapse=''), sep='\t', header=F)
#        colnames(tmp.df) <- c('POS', 'COVERAGE', 'EDITREADS', 'EDITLEVEL', 'REGION')
#        tmp.df <- tmp.df[tmp.df$COVERAGE > 0, ]
#        tmp.df$treatment <- rep(TREATMENT, nrow(tmp.df))
#        treat.df <- rbind(treat.df, tmp.df)
#    }
#    treat.out.df <- aggregate(treat.df[, c('COVERAGE', 'EDITREADS')], by=list(length_df$POS), sum)
#    colnames(treat.out.df) <- c('ID', 'COVERAGE', 'EDITREADS')
#    colnames(treat.out.df) <- c('ID', 'COVERAGE', 'EDITREADS')
#    treat.out.df$EDITLEVEL <- treat.out.df$EDITREADS/treat.out.df$COVERAGE
#    treat.out.df$TREATMENT <- rep(TREATMENT, nrow(treat.out.df))
#    treatment.df <- rbind(treatment.df, treat.out.df)
#}
#
#treatment.df <- treatment.df[treatment.df$COVERAGE > coverage, ]
#filter.site <- tidyr::pivot_wider(treatment.df[, c('ID', 'EDITLEVEL', 'TREATMENT')], values_from='EDITLEVEL', names_from='TREATMENT')
#filter.site$editlevel_pass <- rowSums(filter.site[, unique(treatment.df$TREATMENT)] >= edit_level_thre, na.rm=T)
#filter.site <- filter.site[filter.site$editlevel_pass >= 1, ] # at least one treatment.dfset with editing level >= edit_level_thre
#filter.site$editlevel_high <- rowSums(filter.site[, unique(treatment.df$TREATMENT)] > 0.97, na.rm=T)
#filter.site <- filter.site[filter.site$editlevel_high == 0, c('ID', unique(treatment.df$TREATMENT))] # no treatment.dfsets has editing level >= 0.97
#treatment.df <- treatment.df[treatment.df$ID %in% filter.site$ID]
#write.table(treatment.df, paste(c('replicates_plot', paste(c(out_prefix,'treatment_level.txt'), collapse='.')), collapse='/'), sep='\t', quote=F, row.names=F)


#datasets_id <- c(); more_num <- c()
#for (idx in 1:length(comparison_list)){
#    cmpr <- comparison_list[[idx]]
#    xlab <- cmpr[1]; ylab <- cmpr[2]
#    x.df <- data_region[data_region$SAMPLE == cmpr[1], ]
#    y.df <- data_region[data_region$SAMPLE == cmpr[2], ]
#    combined_ID <- intersect(x.df$ID, y.df$ID)
#    x.df <- x.df[x.df$ID %in% combined_ID, ]; y.df <- y.df[y.df$ID %in% combined_ID, ]
#    x.df$non_EDITREADS <- x.df$COVERAGE - x.df$EDITREADS
#    y.df$non_EDITREADS <- y.df$COVERAGE - y.df$EDITREADS
#    tmp.df <- merge(x.df, y.df, by='ID')
#    tmp.df.x <- tmp.df[, c('EDITREADS.x', 'non_EDITREADS.x', 'EDITREADS.y', 'non_EDITREADS.y')]
#    tmp.df.x$pvalue <- apply(tmp.df.x, 1, function(x) if(any(x < 5)) {fisher.test(matrix(x, nrow=2))$p.value} else{chisq.test(matrix(x,nrow=2))$p.value})
#    tmp.df.x[is.na(tmp.df.x)] <- 1
#    tmp.df.x$ID <- tmp.df$ID
#    tmp.df.x$level <- rep('>0.05', nrow(tmp.df.x))
#    tmp.df.x[tmp.df.x$pvalue <= 0.05 & tmp.df.x$pvalue > 0.01, 'level'] <- '(0.01, 0.05]'
#    tmp.df.x[tmp.df.x$pvalue <= 0.01 & tmp.df.x$pvalue > 0.001, 'level'] <- '(0.001, 0.01]'
#    tmp.df.x[tmp.df.x$pvalue <= 0.001, 'level'] <- '(0, 0.001]'
#    tmp.df.x$EDITLEVEL.x <- x.df[match(tmp.df.x$ID, x.df$ID), 'EDITLEVEL']
#    tmp.df.x$EDITLEVEL.y <- y.df[match(tmp.df.x$ID, y.df$ID), 'EDITLEVEL']
#    wil.pvalue <- round(wilcox.test(tmp.df.x$EDITLEVEL.x, tmp.df.x$EDITLEVEL.y, paired=T, alternative="two.sided")$p.value, 4)
#    r <- round(cor.test(tmp.df.x$EDITLEVEL.x, tmp.df.x$EDITLEVEL.y)$estimate, 4)
#    cor_p <- round(cor.test(tmp.df.x$EDITLEVEL.x, tmp.df.x$EDITLEVEL.y)$p.value, 4)
#    num_x <- round(length(tmp.df.x[tmp.df.x$EDITLEVEL.x > tmp.df.x$EDITLEVEL.y & tmp.df.x$pvalue < 0.05, 'ID']), 4)
#    #datasets_id <- c(datasets_id, cmpr[1]); more_num <- c(more_num, num_x)
#    num_y <- length(tmp.df.x[tmp.df.x$EDITLEVEL.y > tmp.df.x$EDITLEVEL.x & tmp.df.x$pvalue < 0.05, 'ID'])
#    #datasets_id <- c(datasets_id, cmpr[2]); more_num <- c(more_num, num_y)
#    tmp.df.x$level <- factor(tmp.df.x$level, levels=c('>0.05', '(0.01, 0.05]', '(0.001, 0.01]', '(0, 0.001]'))
#    p.tmp.compare <- ggplot(tmp.df.x, aes(x=EDITLEVEL.x, y=EDITLEVEL.y, color=level)) + geom_point(size=1) +
#        geom_abline(intercept=0, slope=1, linetype='dashed', linewidth=2, color='black') +
#        scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1)) +
#        annotate('text', label=paste(c('r=', r, '; p=', cor_p), collapse=''), x=0.8, y=0.2, color='red', size=8) +
#        annotate('text', label=paste(c('n=', num_x), collapse=''), x=0.8, y=0.1, color='red', size=10) +
#        annotate('text', label=paste(c('n=', num_y), collapse=''), x=0.1, y=0.95, color='red', size=10) +
#        annotate('text', label=paste(c('p=', wil.pvalue), collapse=''), x=0.1, y=0.85, color='red', size=10) +
#        scale_color_manual(values=c('#C0C0C0', '#1F78B4', '#FB9A99', '#E31A1C')) + egg::theme_article() +
#        labs(x=xlab, y=ylab, title=paste(c(paste(c(cmpr[1], cmpr[2]), collapse=' vs '), ' n=', length(combined_ID)), collapse='')) +
#        theme(plot.title=element_text(size=15), axis.text.x = element_text(size=16, vjust=0.5), axis.title.x = element_text(size=25), axis.text.y = element_text(size=20), axis.title.y = element_text(size=25), axis.ticks.length.y=unit(0.3, 'cm'), axis.ticks.length.x=unit(0.3, 'cm')) +
#        theme(legend.text=element_text(size=15), legend.title=element_blank(), legend.position='top')
#    if (idx == 1){p.compare <- p.tmp.compare
#    } else {p.compare <- p.compare + p.tmp.compare}
#}
#
#p.compare <- p.compare + plot_layout(ncol=5)
#ggsave(paste(c('replicates_plot', paste(c(out_prefix,'chisq.compare_separate.pdf'), collapse='.')), collapse='/'), p.compare, width=40, height=8*ceiling(length(comparison_list)/5))
