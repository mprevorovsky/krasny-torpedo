##Martin Převorovský prevorov@natur.cuni.cz www.GenoMik.cz

library(RColorBrewer)
chip <- read.delim('./coverage_ChIP-seq/ChIP-seq_gene_coverage.txt', 
                   sep = '\t', header = TRUE, stringsAsFactors = FALSE)
rna <- read.delim('./coverage_RNA-seq/RNA-seq_gene_coverage.txt', 
                   sep = '\t', header = TRUE, stringsAsFactors = FALSE)
which(chip$Name != rna$Name)

pdf('./images/Figure1.pdf')
par(mar = c(5, 5, 5, 5), lwd = 2, cex = 1.5)
options(scipen = -1)
max.chip <- max(chip[, c('cov_wt', 'cov_ko')], na.rm = TRUE)
min.chip <- 10
max.rna <- max(rna[, c('cov_wt', 'cov_ko')], na.rm = TRUE)
min.rna <- 0.1

rbPal <- colorRampPalette(c('yellow', 'red'))

plot(chip$cov_ko, rna$cov_ko,
     main = '',
     xlab = '',
     ylab = '',
     log = 'xy', xlim = c(min.chip, max.chip), ylim = c(min.rna, max.rna),
     type = 'n')
abline(v = c(10, 100, 1000, 10000, 100000), h = c(0.1, 10, 1000, 100000), col = 'gray')
par(new = TRUE)
plot(chip$cov_ko, rna$cov_ko,
     main = 'rnjA KO\nmean of 3 replicates\n(coloured by expression in WT)',
     xlab = 'RNAP occupancy at genes [norm. cov.]',
     ylab = 'gene RNA abundance [norm. cov.]',
     log = 'xy', xlim = c(min.chip, max.chip), ylim = c(min.rna, max.rna), pch = 20,
     col = rbPal(100)[as.numeric(cut(log2(rna$cov_wt + 0.01), breaks = 100))])

plot(chip$cov_wt, rna$cov_wt,
     main = '',
     xlab = '',
     ylab = '',
     log = 'xy', xlim = c(min.chip, max.chip), ylim = c(min.rna, max.rna),
     type = 'n')
abline(v = c(10, 100, 1000, 10000, 100000), h = c(0.1, 10, 1000, 100000), col = 'gray')
par(new = TRUE)
plot(chip$cov_wt, rna$cov_wt,
     main = 'WT\nmean of 3 replicates\n(coloured by expression in WT)',
     xlab = 'RNAP occupancy at genes [norm. cov.]',
     ylab = 'gene RNA abundance [norm. cov.]',
     log = 'xy', xlim = c(min.chip, max.chip), ylim = c(min.rna, max.rna), pch = 20,
     col = rbPal(100)[as.numeric(cut(log2(rna$cov_wt + 0.01), breaks = 100))])
dev.off()