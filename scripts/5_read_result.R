rm(list = ls())
options(stringsAsFactors = FALSE)

resultname <- 'rrbs_clean_data/lmer.RDS'

result <- readRDS(resultname)

range_df <- readRDS('rrbs_clean_data/predictedMeth_range_quality.RDS')
range_df$unit <- rownames(range_df)

library(plyr)
result <- join(result, range_df, by = 'unit', type = "left")
#change var
result$padj <- p.adjust(result$`eGFR:Pr(>|W|)`, method = 'fdr')
result$seqnames <- as.character(result$seqnames)
result <- result[(result$seqnames %in% paste0('chr', 1:22)), ]
saveRDS(result, "rrbs_clean_data/result_range.RDS")


library(openxlsx)
sigresult <- result[result$seqnames %in% paste0('chr', 1:22), ]
#????????X
sigresult <- sigresult[, c('unit', 'seqnames', 'start', 'end', 'eGFR:Pr(>|W|)', 'eGFR:Estimate', 'padj')]
#this place i change the cutoff value to get the total Pvalue to do the HDMT, original cutoff is 0.05 
sigresult <- sigresult[sigresult$'eGFR:Pr(>|W|)' < 1, ]
#????????X
sigresult <- sigresult[order(sigresult$'eGFR:Pr(>|W|)'), ]
#????????X
sigresult <- sigresult[, c('unit', 'seqnames', 'start', 'eGFR:Pr(>|W|)', 'eGFR:Estimate', 'padj')]
colnames(sigresult) <- c('id_row_number', 'Chromosome', 'Position', 'pvalue', 'effectsize', 'padj')
saveRDS(sigresult, file = 'rrbs_clean_data/sigresult.RDS')
write.xlsx(sigresult, file = 'rrbs_clean_data/sigresult.xlsx')

p_cut_v <- c(0.001, 0.01, 0.05, 0.1, 1)

library(openxlsx)
#39??43??49??????????X
for(p_cut in p_cut_v) {
  sigfilename <- gsub('.RDS', paste0('_', p_cut, '.xlsx'), resultname)
  sigfile <- result[result$`eGFR:Pr(>|W|)` < p_cut, ]
  write.xlsx(sigfile, sigfilename)
  
  bedfilename <- gsub('.RDS', paste0('_', p_cut, '.BED'), resultname)
  bedfile <- result[result$`eGFR:Pr(>|W|)` < p_cut, c('seqnames', 'start', 'end', 'unit')]
  write.table(bedfile, bedfilename, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

dir.create('combined-pvalues/data', recursive = TRUE)
bedfilename <- 'combined-pvalues/data/pvals.bed'
bedfile <- result[(result$seqnames %in% paste0('chr', 1:22)), c('seqnames', 'start', 'end', 'eGFR:Pr(>|W|)')]
bedfile$end <- bedfile$start + 1
bedfile$seqnames <- as.character(bedfile$seqnames)
# bedfile$seqnames <- factor(bedfile$seqnames, paste0('chr', 1:22))
bedfile <- bedfile[order(bedfile$seqnames, bedfile$start), ]
colnames(bedfile) <- c('chrom', 'start', 'end', 'p')
bedfile$p <- format(bedfile$p, digits=16, scientific=F)
write.table(bedfile, bedfilename, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


#the following code dont run
library(CMplot)
#62??????????X
cmplotdata <- result[result$seqnames %in% paste0('chr', 1:22), ]
cmplotdata <- cmplotdata[, c('unit', 'seqnames', 'start', 'eGFR:Pr(>|W|)')]
colnames(cmplotdata) <- c('SNP', 'Chromosome', 'Position', 'p-value')
cmplotdata$SNP <- as.factor(cmplotdata$SNP)
cmplotdata$Chromosome <- as.character(cmplotdata$Chromosome)
cmplotdata$Chromosome <- gsub('chr', '', cmplotdata$Chromosome)
cmplotdata$Chromosome <- as.factor(as.numeric(cmplotdata$Chromosome))

cmplotdata <- cmplotdata[cmplotdata$`p-value` > 0 & cmplotdata$`p-value` <= 1, ]

# If the package is version 4.5.1 or 4.5.0, cir.legend is changed to multraits, cir.legend.col is changed to highlight.text.col, and memo is changed to main
CMplot(cmplotdata,plot.type="c",chr.labels=paste("chr",c(1:22),sep=""),r=0.4,multraits=FALSE,
       outward=TRUE,highlight.text.col="black",cir.chr.h=0.3,chr.den.col="black",file="pdf", file.output = TRUE,
       main="",dpi=600,verbose=TRUE)

# If the package is version 4.5.1, Circular-Manhattan.p-value is changed to Cir_Manhtn.p-value
file.rename('Cir_Manhtn.p-value.pdf', 'rrbs_clean_data/Cir_Manhtn.p-value.pdf')

# If the package is version 4.5.1, memo is changed to main
CMplot(cmplotdata,plot.type="q",threshold=1e-6,
       signal.pch=19,signal.cex=1,box=FALSE,multracks=
         FALSE,main="",dpi=600,file = "pdf",file.output=TRUE,verbose=TRUE)

file.rename('QQplot.p-value.pdf', 'rrbs_clean_data/QQplot_p-value.pdf')

