chip <- read.delim('./coverage_ChIP-seq/ChIP-seq_gene_coverage.txt', 
                   sep = '\t', header = TRUE, stringsAsFactors = FALSE)
chip <- cbind(chip, chip$cov_ko / chip$cov_wt)
colnames(chip)[8] <- 'cov_ratio'
rna <- read.csv('./DESeq2results_rnjAKO_vs_WT.txt', 
                header = TRUE, stringsAsFactors = FALSE, row.names = 1)
rna <- rna[!is.na(rna$log2FoldChange), ]
chip <- chip[chip$Name %in% rownames(rna), ]
rna <- rna[rownames(rna) %in% chip$Name, ]
chip <- chip[order(chip$Name), ]
rna <- rna[order(rownames(rna)), ]
which(chip$Name != rownames(rna))

padj.threshold <- 0.05
# class I: >=120% RNAP occupancy (ChIP seq) and sig. upregulated (RNA-seq)
class_I <- chip[which(chip$cov_ratio >= 1.2 & 
                        rna$padj <= padj.threshold & rna$log2FoldChange > 0),
                c('Name', 'locus_tag')]
# class II: 0-80% RNAP occupancy (ChIP seq) and sig. downregulated (RNA-seq)
class_II <- chip[which(chip$cov_ratio <= 0.8 & 
                        rna$padj <= padj.threshold & rna$log2FoldChange < 0),
                c('Name', 'locus_tag')]
# class III: <120% RNAP occupancy in WT (ChIP seq) and sig. upregulated (RNA-seq)
class_III <- chip[which(chip$cov_ratio < 1.2 & 
                        rna$padj <= padj.threshold & rna$log2FoldChange > 0),
                c('Name', 'locus_tag')]
# class IV: >=120% RNAP occupancy (ChIP seq) and sig. downregulated/no change (RNA-seq)
class_IV <- chip[which(chip$cov_ratio >= 1.2 & 
                        ((rna$padj <= padj.threshold & rna$log2FoldChange < 0) |
                          rna$padj > padj.threshold)),
                c('Name', 'locus_tag')]

pdf('./images/classes.pdf')
x.lim <- range(chip$cov_ratio, na.rm = TRUE, finite = TRUE)
y.lim <- range(rna$log2FoldChange, na.rm = TRUE, finite = TRUE)
plot(chip[chip$Name %in% class_I$Name, 'cov_ratio'], 
     rna[rownames(rna) %in% class_I$Name, 'log2FoldChange'],
     log = 'x', xlim = x.lim, ylim = y.lim, 
     xlab = 'RNAP occupancy', ylab = 'mRNA expression')
points(chip[chip$Name %in% class_II$Name, 'cov_ratio'], 
       rna[rownames(rna) %in% class_II$Name, 'log2FoldChange'], 
       col = 'red')
points(chip[chip$Name %in% class_III$Name, 'cov_ratio'], 
       rna[rownames(rna) %in% class_III$Name, 'log2FoldChange'], 
       col = 'blue')
points(chip[chip$Name %in% class_IV$Name, 'cov_ratio'], 
       rna[rownames(rna) %in% class_IV$Name, 'log2FoldChange'], 
       col = 'green')
dev.off()

write.table(class_I, file = 'class_I.txt', quote = FALSE, sep = '\t', row.names = FALSE)
write.table(class_II, file = 'class_II.txt', quote = FALSE, sep = '\t', row.names = FALSE)
write.table(class_III, file = 'class_III.txt', quote = FALSE, sep = '\t', row.names = FALSE)
write.table(class_IV, file = 'class_IV.txt', quote = FALSE, sep = '\t', row.names = FALSE)