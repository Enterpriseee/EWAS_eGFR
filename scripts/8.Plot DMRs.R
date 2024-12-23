options(stringsAsFactors = FALSE)

#load data
result <- read.table("combined-pvalues/data/regions.sig.bed")

#select data
sigresult <- result[result$V7 < 0.05, ]
sigresult <- sigresult[order(sigresult$V7), ]
sigresult <- subset(sigresult, V5>10)
numberOfSignificance <- nrow(sigresult)

#save data
library(xtable)
sigresult_tex <- xtable(sigresult)
print(sigresult_tex, file = "rrbs_clean_data/sigDMR.tex")
write.csv(sigresult, file = "rrbs_clean_data/sigDMR.csv", row.names = FALSE)

#load data
resultname <- "rrbs_clean_data/result_range.RDS"
result <- readRDS(resultname)
cpgresult <- result

#Generate regression plots for each significant region
library(ggplot2)

for(i in 1:nrow(sigresult)) {
  pdf(file = paste0('rrbs_clean_data/DMR_', i, '.pdf'), width = 8.3, height = 4)
  plotdata1 <- sigresult[i, ]
  plotdata2 <- cpgresult[(cpgresult$seqnames == plotdata1$V1) & (cpgresult$start >= plotdata1$V2) & (cpgresult$end <= plotdata1$V3), ]
  print(ggplot(data = plotdata2, mapping = aes(y = plotdata2$`eGFR:Estimate`, x = plotdata2$start)) + geom_smooth() + geom_point() + xlab('BP') + ylab('Coefficient') + ggtitle(paste0('DMR ', plotdata1$V1, ':', plotdata1$V2, '-', plotdata1$V3)) + theme_bw())
  dev.off()
}

library(grid)
library(gridExtra)
fmt_dcimals <- function(decimals=0){
  function(x) sprintf(paste0("%-", decimals, "s"), x)
}

#Display all significant area maps on multiple pages
pdf(file = paste0("rrbs_clean_data/DMRs.pdf"), width = 28, height = 16.5)
col_per_page <- 4
row_per_page <- 4
num_per_page <- col_per_page * row_per_page
num_page <- ceiling(nrow(sigresult)/num_per_page)
for(npage in 1:num_page) {
  startplot <- (npage -1) * num_per_page + 1
  if((startplot + num_per_page) > nrow(sigresult)) {
    endplot <- nrow(sigresult)
  } else {
    endplot <- startplot + num_per_page - 1
  }
  cat(paste0("start: ", startplot, ", end: ", endplot, "\n"))
  plot_list <- list()
  for(i in startplot:endplot) {
    plotdata1 <- sigresult[i, ]
    plotdata2 <- cpgresult[(cpgresult$seqnames == plotdata1$V1) & (cpgresult$start >= plotdata1$V2) & (cpgresult$end <= plotdata1$V3), ]
    plot_list[[i-startplot+1]] <- ggplot(data = plotdata2, mapping = aes(y = `eGFR:Estimate`, x = start)) + geom_smooth() + geom_point() + xlab('BP') + ylab('Coefficient') + ggtitle(bquote(bold(.(LETTERS[i-startplot+1])) * ' DMR ' * .(plotdata1$V1) * ':' * .(plotdata1$V2) * '-' * .(plotdata1$V3))) + theme_bw() + scale_y_continuous(labels = fmt_dcimals(4)) + theme(axis.text.x = element_text(angle = 0, hjust = 1))
  }
  do.call('grid.arrange', c(plot_list, ncol = col_per_page, nrow = row_per_page))
}
dev.off()

library(grid)
library(gridExtra)
fmt_dcimals <- function(decimals=0){
  function(x) sprintf(paste0("%-", decimals, "s"), x)
}

#Another format for multi-page presentation
pdf(file = paste0("rrbs_clean_data/DMRsr.pdf"), width = 15, height = 11.7)
col_per_page <- 3
row_per_page <- 4
num_per_page <- col_per_page * row_per_page
num_page <- ceiling(nrow(sigresult)/num_per_page)
for(npage in 1:num_page) {
  startplot <- (npage -1) * num_per_page + 1
  if((startplot + num_per_page) > nrow(sigresult)) {
    endplot <- nrow(sigresult)
  } else {
    endplot <- startplot + num_per_page - 1
  }
  cat(paste0("start: ", startplot, ", end: ", endplot, "\n"))
  plot_list <- list()
  for(i in startplot:endplot) {
    plotdata1 <- sigresult[i, ]
    plotdata2 <- cpgresult[(cpgresult$seqnames == plotdata1$V1) & (cpgresult$start >= plotdata1$V2) & (cpgresult$end <= plotdata1$V3), ]
    plot_list[[i-startplot+1]] <- ggplot(data = plotdata2, mapping = aes(y = `eGFR:Estimate`, x = start)) + geom_smooth() + geom_point() + xlab('BP') + ylab('Coefficient') + ggtitle(bquote(bold(.(LETTERS[i-startplot+1])) * ' DMR ' * .(plotdata1$V1) * ':' * .(plotdata1$V2) * '-' * .(plotdata1$V3))) + theme_bw() + scale_y_continuous(labels = fmt_dcimals(4)) + theme(axis.text.x = element_text(angle = , hjust = 1))
  }
  do.call('grid.arrange', c(plot_list, ncol = col_per_page, nrow = row_per_page))
}
dev.off()

