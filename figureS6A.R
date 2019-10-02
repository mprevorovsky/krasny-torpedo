# https://stat.ethz.ch/pipermail/bioconductor/2008-October/024669.html
getAttributeField <- function (x, field, attrsep = ";") {
  s = strsplit(x, split = attrsep, fixed = TRUE)
  sapply(s, function(atts) {
    a = strsplit(atts, split = "=", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv <- a[[m]][2]
    }
    else {
      rv <- as.character(NA)
    }
    return(rv)
  })
}

gffRead <- function(gffFile, nrows = -1) {
  cat("Reading ", gffFile, ": ", sep="")
  gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
                   header=FALSE, comment.char="#", nrows = nrows,
                   colClasses=c("character", "character", "character", "integer",  
                                "integer",
                                "character", "character", "character", "character"))
  colnames(gff) = c("seqname", "source", "feature", "start", "end",
                    "score", "strand", "frame", "attributes")
  cat("found", nrow(gff), "rows with classes:",
      paste(sapply(gff, class), collapse=", "), "\n")
  stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
  return(gff)
}
##

##Martin Převorovský prevorov@natur.cuni.cz www.GenoMik.cz
genelist <- c('ldh', 'ydbL', 'ymcB', 'bsrB', 'veg', 'rpmI', 'hag', 'purR', 'pyrR')

gff_file <- './genome/GCF_000009045.1_ASM904v1_genomic_reverse_SYS.gff'
gff <- gffRead(gff_file)
gff <- gff[gff$feature == 'gene', ]
gff$Name <- getAttributeField(gff$attributes, 'Name')
gff$locus_tag <- getAttributeField(gff$attributes, 'locus_tag')
gff[which(gff$locus_tag == 'BSU_misc_RNA_32'), 'locus_tag'] <- 'bsrB'

ko <- read.csv('./coverage_ChIP-seq/KO.input-normalized.mean.bedgraph', head = FALSE, stringsAsFactors = FALSE, skip = 1, sep = '\t')
wt <- read.csv('./coverage_ChIP-seq/wt.input-normalized.mean.bedgraph', head = FALSE, stringsAsFactors = FALSE, skip = 1, sep = '\t')
data <- rbind(wt, ko)

for (i in genelist) {
  start <- gff[which(gff$locus_tag == i) - 1, 'start']
  end <- gff[which(gff$locus_tag == i) + 1, 'end']
  chr <- gff[which(gff$locus_tag == i), 'seqname']
  ymax <- max(data[data[, 2] >= start & data[, 2] <= end, 4], na.rm = TRUE)
  ymin <- (-1 * ymax) / 4

  pdf(file = paste0('./images/', start, '-', end, '_', i, '.pdf'), title = paste0(i, '.pdf'))
  par(mfrow = c(2, 1), mar = c(3, 4, 4, 2))
  
  plot(wt[which(wt[, 1] == chr & wt[, 2] >= start & wt[, 2] <= end), c(2, 4)],
       type = 'h', col = 'blue', ylim = c(ymin, ymax), 
       main = 'WT', xlab = '', ylab = 'normalized coverage', cex.axis = 1.2, cex.lab = 1.2)
  features = gff[gff$seqname == chr & (gff$start <= end & gff$end >= start), ]
  for (i in 1:nrow(features)){
    if (features[i, 'strand'] == '-'){
      arrows(features[i, 'start'], ymin, features[i, 'end'], ymin, lwd = 4, length = 0.1)
    }
    else{
      arrows(features[i, 'end'], ymin, features[i, 'start'], ymin, lwd = 4, length = 0.03)
    }
  }
  text(labels = features$locus_tag, x = features$start + (features$end - features$start) / 2, 
       y = ymin + ymax / 10, cex = 1.5)
  
  plot(ko[which(ko[, 1] == chr & ko[, 2] >= start & ko[, 2] <= end), c(2, 4)],
       type = 'h', col = 'red', ylim = c(ymin, ymax), 
       main = 'rnjA', xlab = '', ylab = 'normalized coverage', cex.axis = 1.2, cex.lab = 1.2)
  features = gff[gff$seqname == chr & (gff$start <= end & gff$end >= start), ]
  for (i in 1:nrow(features)){
    if (features[i, 'strand'] == '-'){
      arrows(features[i, 'start'], ymin, features[i, 'end'], ymin, lwd = 4, length = 0.1)
    }
    else{
      arrows(features[i, 'end'], ymin, features[i, 'start'], ymin, lwd = 4, length = 0.03)
    }
  }
  text(labels = features$locus_tag, x = features$start + (features$end - features$start) / 2, 
       y = ymin + ymax / 10, cex = 1.5)
  
  dev.off()  
}