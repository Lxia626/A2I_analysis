library(ggplot2)

args <- commandArgs(TRUE)
FUNCTION <- args[1]

if (FUNCTION == "percent_single") {
    data_prefix <- args[2]
    out_prefix <- args[3]
    # mutation type of single sample 
    data <- read.table(paste(data_prefix, "edit_type.txt", sep="."), header=TRUE, sep="\t")
    sample_id <- strsplit(data_prefix, "\\.")[[1]][1]
    data$TYPE <- gsub("A->T", "others", data$TYPE)
    data$TYPE <- gsub("A->C", "others", data$TYPE)
    data$TYPE <- gsub("C->G", "others", data$TYPE)
    data$TYPE <- gsub("C->A", "others", data$TYPE)
    data$TYPE <- gsub("G->T", "others", data$TYPE)
    data$TYPE <- gsub("G->C", "others", data$TYPE)
    data$TYPE <- gsub("T->A", "others", data$TYPE)
    data$TYPE <- gsub("T->G", "others", data$TYPE)
    data$TYPE <- factor(data$TYPE, levels=c("A->G", "T->C", "C->T", "G->A", "others"))
    data <- data[match(levels(data$TYPE), data$TYPE), ]
    pdf(paste(out_prefix, "edit_type", "pdf", sep="."), width=6, height=6)
    pie(data$PERCENT, labels=with(data, paste0(TYPE, " (", round(PERCENT, 4) * 100, "%)")), main=paste(sample_id, "editing_type", sep="."))
    dev.off()

    # motif (up: G-enriched, down: T/U-enriched) 
    motif.data <- read.table(paste(data_prefix, "motif.txt", sep="."), header=TRUE, sep="\t")
    motif.data <- motif.data[, c("UP", "DOWN")]
    up.data <- aggregate(motif.data$UP, by=list(motif.data$UP), length)
    colnames(up.data) <- c("BASE", "COUNT")
    down.data <- aggregate(motif.data$DOWN, by=list(motif.data$DOWN), length)
    colnames(down.data) <- c("BASE", "COUNT")
    motif <- dplyr::left_join(up.data, down.data, by="BASE")
    colnames(motif) <- c("BASE", "UP", "DOWN")
    motif.df <- tidyr::pivot_longer(motif, cols=2:3, names_to="POS", values_to="COUNT")
    motif.df$BASE <- factor(motif.df$BASE, levels=c("A", "C", "G", "T"))
    motif.df$POS <- factor(motif.df$POS, levels=c("UP", "DOWN"))
    p.motif <- ggplot(motif.df, aes(x=POS, y=COUNT, fill=BASE)) + geom_bar(stat="identity") + egg::theme_article() +
        theme(axis.text.x = element_text(size=15), axis.title.x = element_text(size=20), axis.text.y = element_text(size=15), axis.title.y = element_text(size=20))
    ggsave(paste(out_prefix, "motif", "pdf", sep="."), width=10, height=6)

} else if (FUNCTION == 'percent_multiple') {

    library(patchwork)
    library(RColorBrewer)
    sample_names <- args[2]
    sample_order <- strsplit(sample_names, ",")[[1]]
    outFile <- args[3]

    type.plot <- data.frame(); labels <- c()
    for (SAMPLE in sample_order){
        data_alu <- read.table(paste(SAMPLE, 'Alu.edit_type.txt', sep='.'), header=TRUE, sep='\t')
        SUM <- sum(data_alu$COUNT)
        A2G_ratio <- data_alu[data_alu$TYPE == 'A->G', 'PERCENT']
        T2C_ratio <- data_alu[data_alu$TYPE == 'T->C', 'PERCENT']
        C2T_ratio <- data_alu[data_alu$TYPE == 'C->T', 'PERCENT']
        G2A_ratio <- data_alu[data_alu$TYPE == 'G->A', 'PERCENT']
        A2G_T2C_ratio <- sum(data_alu[data_alu$TYPE %in% c('A->G', 'T->C'), 'PERCENT'])
        label <- paste(SUM, expression('\n'), A2G_ratio, expression('\n'), T2C_ratio, expression('\n'), C2T_ratio, expression('\n'), G2A_ratio, expression('\n'), A2G_T2C_ratio)
        labels <- c(labels, label)
        data_alu_other <- data.frame(TYPE='other', COUNT=sum(data_alu[!data_alu$TYPE %in% c('A->G', 'T->C', 'C->T', 'G->A'), ]$COUNT), PERCENT=sum(data_alu[!data_alu$TYPE %in% c('A->G', 'T->C', 'C->T', 'G->A'), ]$PERCENT))
        data_alu <- rbind(data_alu, data_alu_other)
        data_alu <- data_alu[data_alu$TYPE %in% c('A->G', 'T->C', 'C->T', 'G->A', 'other'), ]
        data_alu$sample <- rep(SAMPLE, length(data_alu$TYPE)); data_alu$rep <- rep('Alu', length(data_alu$TYPE))
        type.plot <- rbind(type.plot, data_alu)
        data_rep <- read.table(paste(SAMPLE, 'NoAlu_rep.edit_type.txt', sep='.'), header=TRUE, sep='\t')
        data_rep$sample <- rep(SAMPLE, length(data_rep$TYPE)); data_rep$rep <- rep('NoAlu_rep', length(data_rep$TYPE))
        type.plot <- rbind(type.plot, data_rep)
        data_other <- read.table(paste(SAMPLE, 'other.edit_type.txt', sep='.'), header=TRUE, sep='\t')
        data_other$sample <- rep(SAMPLE, length(data_other$TYPE)); data_other$rep <- rep('other', length(data_other$TYPE))
        type.plot <- rbind(type.plot, data_other)
    }
    count_text <- data.frame(sample=strsplit(sample_names, ",")[[1]], LABEL=labels, Y=rep(0.8, length(labels)))
    alu.plot <- type.plot[type.plot$rep == 'Alu', ]
    alu.plot$TYPE <- factor(alu.plot$TYPE, levels=c('other', 'G->A', 'C->T', 'T->C', 'A->G'))
    rep.plot <- type.plot[type.plot$TYPE == 'A->G', ]
    rep.plot.sum <- aggregate(rep.plot$COUNT, by=list(rep.plot$sample, rep.plot$TYPE), sum)
    colnames(rep.plot.sum) <- c('sample', 'TYPE', 'Sum')
    rep.plot <- dplyr::left_join(rep.plot, rep.plot.sum, by='sample')
    rep.plot$PERCENT <- rep.plot$COUNT/rep.plot$Sum
    rep.plot$rep <- factor(rep.plot$rep, levels=c('other', 'NoAlu_rep', 'Alu'))
    
    alu.plot$sample <- factor(alu.plot$sample, levels=sample_order)
    p.type <- ggplot() + geom_bar(alu.plot, mapping=aes(x=sample, y=PERCENT, fill=TYPE), stat='identity', width=0.8) + scale_y_continuous(expand = expansion(mult = c(0,0))) +
        egg::theme_article() + scale_fill_manual(values=c('#D3D3D3', '#87CEFA', '#6495ED', '#FFDAB9', '#FF7F50')) + 
        theme(axis.text.x=element_text(size=15, angle=90), axis.text.y=element_text(size=15), axis.title.x=element_text(size=18), axis.title.y=element_text(size=18), legend.title=element_blank(), legend.text=element_text(size=12)) +
        geom_text(count_text, mapping=aes(x=sample, y=Y, label=LABEL), size=4)
    rep.plot$sample <- factor(rep.plot$sample, levels=sample_order)
    p.rep <- ggplot(rep.plot, aes(x=sample, y=PERCENT)) + geom_bar(aes(fill=rep), stat='identity', width=0.8) + scale_y_continuous(expand = expansion(mult = c(0,0))) +
        egg::theme_article() + scale_fill_manual(values=c('#8DD3C7', '#FFFFB3', '#BEBADA')) + 
        theme(axis.text.x=element_text(size=15, angle=90), axis.text.y=element_text(size=15), axis.title.x=element_text(size=18), axis.title.y=element_text(size=18), legend.title=element_blank(), legend.text=element_text(size=12))
    p.type_rep <- p.type + p.rep
    ggsave(outFile, p.type_rep, width=25, height=8)

}
