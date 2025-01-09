#! /usr/bin/env Rscript

.libPaths(c('/public/home/Songlab/xial/miniconda3/lib/R/library'))
# .libPaths()
library(qpdf)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
level_dir <- args[1] # e.g tmp_dir/, note that '/' is needed.
coverage <- as.numeric(args[2])
na_thresh <- as.numeric(args[3])
editlevel <- args[4]
out_prefix <- args[5]
min_editlevel <- as.numeric(strsplit(editlevel, ',')[[1]][1])
max_editlevel <- as.numeric(strsplit(editlevel, ',')[[1]][2])
list_samples <- list.files(level_dir)
sample_list <- sapply(list_samples, function(x) strsplit(x,'\\.')[[1]][1])
dir.create('replicates_plot')

level_df <- read.table(paste(c(level_dir, list_samples[1]), collapse=''), header=F, sep='\t')
site_region <- level_df[, c('V1', 'V5')]
colnames(site_region) <- c('ID', 'region')
level_df$V4 <- gsub('N/A', NA, level_df$V4)
sample_id <- strsplit(list_samples[1], '\\.')[[1]][1]
level_df[level_df$V2 < coverage, 'V4'] <- NA
level_df <- level_df[, c('V1', 'V4')]
level_df$V4 <- as.numeric(level_df$V4)
colnames(level_df) <- c('ID', sample_id)
print(paste(c(list_samples[1], dim(level_df)[1]), collapse=': '))
level_df <- unique(level_df)

for (FILE in list_samples[2:length(list_samples)]){
    sample_id <- strsplit(FILE, '\\.')[[1]][1]
    tmp.df <- read.table(paste(c(level_dir,FILE), collapse=''), header=F, sep='\t')
    tmp.df$V4 <- gsub('N/A', NA, tmp.df$V4)
    tmp.df[tmp.df$V2 < coverage, 'V4'] <- NA
    tmp.df <- tmp.df[, c('V1', 'V4')]
    tmp.df$V4 <- as.numeric(tmp.df$V4)
    colnames(tmp.df) <- c('ID', sample_id)
    print(paste(c(list_samples[1], dim(level_df)[1]), collapse=': '))
    tmp.df <- unique(tmp.df)
    level_df <- left_join(level_df, tmp.df, by='ID')
}

row.names(level_df) <- level_df[, 'ID']
level_df.plot <- level_df[, sample_list]
if (min_editlevel == 0){
    level_df.plot$editlevel_pass <- rowSums(level_df.plot == 0, na.rm=T)
    level_df.plot <- level_df.plot[level_df.plot$editlevel_pass < length(sample_list), sample_list] # omit sites where all samples do not edit.
} else {
    level_df.plot$editlevel_pass <- rowSums(level_df.plot > min_editlevel, na.rm=T)
    if (na_thresh >= 1) {level_df.plot <- level_df.plot[level_df.plot$editlevel_pass >= na_thresh, sample_list] # at least one sample shows editing level >= 1%, according to Gabay et al ., Nat Commun, 2022
    } else {level_df.plot <- level_df.plot[level_df.plot$editlevel_pass/length(sample_list) >= na_thresh, sample_list]} # at least na_thresh percent samples show editing level >= 1%
}
level_df.plot$editlevel_high <- rowSums(level_df.plot > max_editlevel, na.rm=T)
#level_df.plot <- level_df.plot[level_df.plot$editlevel_high < length(sample_list), sample_list]
level_df.plot <- level_df.plot[level_df.plot$editlevel_high == 0, sample_list]
#level_df.plot$max_level <- apply(level_df.plot[, sample_list], 1, function(x) max(x, na.rm=T)) # discard those with editing level >= 97%, according to Song et al., Genom Res, 2018, 2022
#level_df.plot <- level_df.plot[level_df.plot$max_level < max_editlevel, sample_list]
#print('aa')
#print(head(level_df[level_df$na_num < length(sample_list), c(sample_list, 'editlevel_pass', 'max_level', 'na_num')]))
#print('bb')
#print(head(level_df[level_df$na_num == length(sample_list), c(sample_list, 'editlevel_pass', 'max_level', 'na_num')]))
#print('cc')
#print(head(level_df[!(is.na(level_df$editlevel_pass)) & level_df$editlevel_pass >= na_thresh, c(sample_list, 'editlevel_pass', 'max_level', 'na_num')]))
#print('dd')
#print(head(level_df[!(is.na(level_df$max_level)) & level_df$max_level < max_editlevel, c(sample_list, 'editlevel_pass', 'max_level', 'na_num')]))
#print('ee')
#print(head(level_df[level_df$na_num < length(sample_list) & level_df$max_level < max_editlevel, c(sample_list, 'editlevel_pass', 'max_level', 'na_num')]))
level_df.plot <- cbind(ID=rownames(level_df.plot), level_df.plot)
row.names(level_df.plot) <- NULL
write.table(level_df.plot, paste(c('replicates_plot', paste(c(out_prefix,'filtered.edit_level_noRegion.txt'),collapse='.')), collapse='/'), quote=F, sep='\t', row.names=F)
level_df.out <- left_join(level_df.plot, site_region, by='ID')
write.table(level_df.out, paste(c('replicates_plot', paste(c(out_prefix,'filtered.edit_level.txt'),collapse='.')), collapse='/'), quote=F, sep='\t', row.names=F)

### analysis of replicates for each datasets
level_df.out <- read.table(paste(c('replicates_plot', paste(c(out_prefix,'filtered.edit_level.txt'),collapse='.')), collapse='/'), quote='', header=T, sep='\t')
## Venn plot for overlapped editing sites of single datasets
# library(RColorBrewer)
library(ggplot2)
theme_set(theme_classic(base_family = "Arial"))
library(GGally)
venn.site <- function(list_data, category_names, prefix, region){
    if (length(category_names) >= 3 & length(category_names) <= 5){
        library(VennDiagram)
        library(grid)
        venn.plot <- venn.diagram(x=list_data, category.names=category_names, fill=brewer.pal(length(category_names),'Pastel1'), print.mode=c("raw","percent"),
            main=paste(c(region,'pdf'), collapse='.'), main.cex=1.5, 
            cex=1.5, cat.cex=2, cat.dist=0.1, cat.fontface='bold', filename=NULL)
        pdf(paste(c(prefix, 'pdf'), collapse='.'), width=8, height=8)
        grid.draw(venn.plot)
        dev.off()
    } else if (length(category_names) < 3) {
        library(VennDiagram)
        library(grid)
        venn.plot <- venn.diagram(x=list_data, category.names=category_names, print.mode=c("raw","percent"),
            main=paste(c(region,'pdf'), collapse='.'), main.cex=1.5, 
            cex=1.5, cat.cex=2, cat.dist=0.1, cat.fontface='bold', filename=NULL)
        pdf(paste(c(prefix, 'pdf'), collapse='.'), width=8, height=8)
        grid.draw(venn.plot)
        dev.off()
    } else {
        .libPaths(c('/public/home/Songlab/xial/miniconda3/envs/DE_analysis/lib/R/library', .libPaths()))
        library(venn)
        pdf(paste(c(prefix, 'pdf'), collapse='.'), width=8, height=8)
        venn(list_data, snames=category_names, sncs=1, ilcs=0.9, zcolor=brewer.pal(length(category_names), 'Set2'), box=F)
        dev.off()
    }
}

# number of overlapped editing sites
#pdf_file_list <- c()
#for (region in c('ALU', 'nonALU_rep', 'other')){
#    tmp.level_df <- level_df.out[level_df.out$region == region, ]
#    site_list <- list()
#    for (x in sample_list){
#        site_list <- append(site_list, list(na.omit(tmp.level_df[, c('ID', x)])[, 'ID']))
#    }
#    venn.site(site_list, sample_list, paste(c('replicates_plot', paste(c(out_prefix,region,'overlapped_site.venn'),collapse='.')), collapse='/'), region)
#    pdf_file_list <- c(pdf_file_list, paste(c('replicates_plot', paste(c(out_prefix,region,'overlapped_site.venn.pdf'),collapse='.')), collapse='/'))
#}
#pdf_combine(pdf_file_list, output = paste(c('replicates_plot', paste(c(out_prefix,'joined.overlapped_site.venn.pdf'),collapse='.')), collapse='/'))
#site_list_total <- list()
#for (x in sample_list){
#    site_list_total <- append(site_list_total, list(na.omit(level_df.plot[, c('ID', x)])[, 'ID']))
#}
#venn.site(site_list_total, sample_list, paste(c('replicates_plot', paste(c(out_prefix,'all.overlapped_site.venn'),collapse='.')), collapse='/'), 'all')

kruskal_calculate <- function(inputData, Y) {
    single_compare <- agricolae::kruskal(inputData$EDIT_LEVEL, inputData$SAMPLE_ID)$groups
    rownamemar <- row.names(single_compare)
    out_kruskal.single <- data.frame(rownamemar, single_compare$groups) %>% filter(rownamemar!="NA")
    colnames(out_kruskal.single) <- c('SAMPLE_ID', 'marker')
    out_kruskal.single$Y <- Y
    return(out_kruskal.single)
}

# editing level boxplot
level_df.plot_cmb <- level_df.out[!(duplicated(level_df.out$ID)), ]
max_Y <- max(apply(level_df.plot_cmb[,2:(ncol(level_df.plot_cmb)-1)], 2, function(x) fivenum(x, na.rm=T)[4]+1.5*(fivenum(x, na.rm=T)[4]-fivenum(x, na.rm=T)[2])))
level_df.box <- tidyr::pivot_longer(level_df.plot_cmb, col=2:8, values_to='EDIT_LEVEL', names_to='SAMPLE_ID')
out_kruskal.single <- kruskal_calculate(level_df.box, max(na.omit(level_df.box$EDIT_LEVEL))*1.05)
out_kruskal.single$SAMPLE_ID <- factor(out_kruskal.single$SAMPLE_ID, levels=sort(unique(out_kruskal.single$SAMPLE_ID)))
level_df.box$SAMPLE_ID <- factor(level_df.box$SAMPLE_ID, levels=sort(unique(level_df.box$SAMPLE_ID)))
level_df.box$region <- factor(level_df.box$region, levels=c('ALU', 'nonALU_rep', 'other'))
p.box_jit <- ggplot(level_df.box, aes(x=SAMPLE_ID, y=EDIT_LEVEL)) + 
    geom_point(size=1, position = position_jitter(width = 0.3, height = 0), alpha=0.25) +
    geom_boxplot(outlier.shape=NA, fill=NA, color='red') +
    # scale_y_continuous(limits=c(0, 1.1)) +
    # scale_x_discrete(labels=c('PE101_1', 'PE101_2', 'PE76_1', 'PE76_2', 'PE76_3', 'PE76_4', 'PE76_5')) +
    geom_text(data=out_kruskal.single, mapping=aes(x=SAMPLE_ID, y=Y, label=marker), size=5, position=position_dodge(0.6)) +
    egg::theme_article() +
    theme(axis.text=element_text(size=10, color='black'), axis.title=element_text(size=15, color='black'), axis.ticks.length=unit(3,'mm'), axis.line=element_line(linewidth=0.25, color='black'))
ggsave(paste(c('replicates_plot', paste(c(out_prefix,'na_level_box_jitter.pdf'),collapse='.')), collapse='/'), width=7, height=5)

out_kruskal.single <- kruskal_calculate(level_df.box, max_Y*1.05)
out_kruskal.single$SAMPLE_ID <- factor(out_kruskal.single$SAMPLE_ID, levels=sort(unique(out_kruskal.single$SAMPLE_ID)))
p.box<- ggplot(level_df.box, aes(x=SAMPLE_ID, y=EDIT_LEVEL)) + 
    geom_boxplot(outlier.shape=NA, fill=NA, color='red') +
    coord_cartesian(ylim = c(0,max_Y*1.05)) +
    # scale_x_discrete(labels=c('PE101_1', 'PE101_2', 'PE76_1', 'PE76_2', 'PE76_3', 'PE76_4', 'PE76_5')) +
    geom_text(data=out_kruskal.single, mapping=aes(x=SAMPLE_ID, y=Y, label=marker), size=5, position=position_dodge(0.6)) +
    egg::theme_article() +
    theme(axis.text=element_text(size=10, color='black'), axis.title=element_text(size=15, color='black'), axis.ticks.length=unit(3,'mm'), axis.line=element_line(linewidth=0.25, color='black'))
ggsave(paste(c('replicates_plot', paste(c(out_prefix,'na_level_box.pdf'),collapse='.')), collapse='/'), width=7, height=5)
q(save='no')
# editing level paired analysis
level_df.plot_cmb <- level_df.out[!(duplicated(level_df.out$ID)), c(sample_list, 'region')]
my_fn <- function(data, mapping, pts=list(), abline=list(), ...){
    ggplot(data = data, mapping = mapping,) +
         do.call(geom_point, pts) +
         do.call(geom_abline, abline)
}
pair_plot <- function(level_df.plot_cmb, region){
    if (region %in% c('ALU', 'nonALU_rep', 'other')) {
        level_df.plot_cmb <- level_df.plot_cmb[level_df.plot_cmb$region == region, ]
    }
    pair_level.p <- ggpairs(na.omit(level_df.plot_cmb), columns=1:length(sample_list), mapping=aes(color=region), upper=list(continuous=wrap('cor', size=5)), lower=list(continuous=wrap(my_fn, pts=list(mapping=aes(color=region), size=1), abline=list(intercept=0, slope=1, linetype='dashed', linewidth=2, color='black'))), cardinality_threshold = NULL) +
        egg::theme_article() + labs(title=paste(c('ggpair.for.', region, ' (n=', dim(na.omit(level_df.plot_cmb))[1], ')'), collapse='')) +
        theme(axis.text.x = element_text(size=15), axis.title.x = element_text(size=20), axis.text.y = element_text(size=15), axis.title.y = element_text(size=20), axis.ticks.length.x=unit(0.3, 'cm'), axis.ticks.length.y=unit(0.3, 'cm')) +
        theme(strip.text.y=element_text(size=20), strip.text.x=element_text(size=20), plot.title=element_text(size=25))
    ggsave(paste(c('replicates_plot', paste(c(out_prefix,region,'no_na_edit_level.pair.pdf'),collapse='.')), collapse='/'), pair_level.p, width=13, height=13, limitsize=F)
}
pair_plot_na <- function(level_df.plot_cmb, region){
    if (region %in% c('ALU', 'nonALU_rep', 'other')) {
        level_df.plot_cmb <- level_df.plot_cmb[level_df.plot_cmb$region == region, ]
    }
    pair_level.p <- ggpairs(level_df.plot_cmb, columns=1:length(sample_list), mapping=aes(color=region), upper=list(continuous=wrap('cor', size=5)), lower=list(continuous=wrap(my_fn, pts=list(mapping=aes(color=region), size=1), abline=list(intercept=0, slope=1, linetype='dashed', linewidth=2, color='black'))), cardinality_threshold = NULL) +
        egg::theme_article() + labs(title=paste(c('ggpair.for.', region, ' (n=', dim(level_df.plot_cmb)[1], ')'), collapse='')) +
        theme(axis.text.x = element_text(size=15), axis.title.x = element_text(size=20), axis.text.y = element_text(size=15), axis.title.y = element_text(size=20), axis.ticks.length.x=unit(0.3, 'cm'), axis.ticks.length.y=unit(0.3, 'cm')) +
        theme(strip.text.y=element_text(size=20), strip.text.x=element_text(size=20), plot.title=element_text(size=25))
    ggsave(paste(c('replicates_plot', paste(c(out_prefix,region,'na_edit_level.pair.pdf'),collapse='.')), collapse='/'), pair_level.p, width=13, height=13, limitsize=F)
}

pair_plot(level_df.plot_cmb, 'ALU')
pair_plot(level_df.plot_cmb, 'nonALU_rep')
pair_plot(level_df.plot_cmb, 'other')
pair_plot(level_df.plot_cmb, 'all')
pair_plot_list <- c(paste(c('replicates_plot', paste(c(out_prefix,'ALU','no_na_edit_level.pair.pdf'),collapse='.')), collapse='/'))
pair_plot_list <- c(pair_plot_list, paste(c('replicates_plot', paste(c(out_prefix,'nonALU_rep','no_na_edit_level.pair.pdf'),collapse='.')), collapse='/'))
pair_plot_list <- c(pair_plot_list, paste(c('replicates_plot', paste(c(out_prefix,'other','no_na_edit_level.pair.pdf'),collapse='.')), collapse='/'))
pdf_combine(pair_plot_list, output = paste(c('replicates_plot', paste(c(out_prefix,'joined.no_na_edit_level.pair.pdf'),collapse='.')), collapse='/'))
pair_plot_na(level_df.plot_cmb, 'ALU')
pair_plot_na(level_df.plot_cmb, 'nonALU_rep')
pair_plot_na(level_df.plot_cmb, 'other')
pair_plot_na(level_df.plot_cmb, 'all')
pair_plot_list_na <- c(paste(c('replicates_plot', paste(c(out_prefix,'ALU','na_edit_level.pair.pdf'),collapse='.')), collapse='/'))
pair_plot_list_na <- c(pair_plot_list_na, paste(c('replicates_plot', paste(c(out_prefix,'nonALU_rep','na_edit_level.pair.pdf'),collapse='.')), collapse='/'))
pair_plot_list_na <- c(pair_plot_list_na, paste(c('replicates_plot', paste(c(out_prefix,'other','na_edit_level.pair.pdf'),collapse='.')), collapse='/'))
pdf_combine(pair_plot_list_na, output = paste(c('replicates_plot', paste(c(out_prefix,'joined.na_edit_level.pair.pdf'),collapse='.')), collapse='/'))
