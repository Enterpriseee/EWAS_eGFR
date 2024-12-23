#Two functions are defined: bsraw_matrix() and bsrel_matrix() to process the RRBS data and generate two files (containing location information and methylation level data, respectively).

options(stringsAsFactors = FALSE)

#define bsraw_matrix()
bsraw_matrix <- function(filename) {
  
  #library software and load data
  require(BiSeq)
  bsraw <- readRDS(filename)
  
  #Preprocess the data
  mgranges <- rowRanges(bsraw)
  mgranges <- as.data.frame(mgranges)
  rownames(mgranges) <- paste0('r', 1:nrow(mgranges))
  
  mgranges_new <- rbind.data.frame(mgranges[2:nrow(mgranges), ], mgranges[1, ])
  idx_keep <- rownames(mgranges)[!((mgranges$seqnames == mgranges_new$seqnames) & (mgranges$start == (mgranges_new$start - 1)))]
  mgranges <- mgranges[idx_keep, ]
  
  mtotal <- totalReads(bsraw)
  mmeth <- methReads(bsraw)
  
  #calculat methylation levels
  mbeta <- mmeth/mtotal
  mbeta[which(is.nan(mbeta), arr.ind = TRUE)] <- NA
  rownames(mbeta) <- paste0('r', 1:nrow(mbeta))
  
  #select CpGs
  mbeta <- mbeta[idx_keep, ]
  
  #save data
  file_range <- gsub('.RDS', '_range.RDS', filename)
  file_beta <- gsub('.RDS', '_beta.RDS', filename)
  
  saveRDS(mgranges, file_range)
  saveRDS(mbeta, file_beta)
}

#run bsraw_matrix
bsraw_matrix('rrbs_clean_data/rrbs.clust.lim.RDS')

#define bsrel_matrix
bsrel_matrix <- function(filename) {
  
  #library software and load data
  require(BiSeq)
  bsrel <- readRDS(filename)
  
  #Preprocess the data
  mgranges <- rowRanges(bsrel)
  mgranges <- as.data.frame(mgranges)
  rownames(mgranges) <- paste0('r', 1:nrow(mgranges))
  
  mgranges_new <- rbind.data.frame(mgranges[2:nrow(mgranges), ], mgranges[1, ])
  idx_keep <- rownames(mgranges)[!((mgranges$seqnames == mgranges_new$seqnames) & (mgranges$start == (mgranges_new$start - 1)))]
  mgranges <- mgranges[idx_keep, ]
  
  #extract methylation levels at each site
  mbeta <- methLevel(bsrel)
  rownames(mbeta) <- paste0('r', 1:nrow(mbeta))
  mbeta <- mbeta[idx_keep, ]
  
  #save data
  file_range <- gsub('.RDS', '_range.RDS', filename)
  file_beta <- gsub('.RDS', '_beta.RDS', filename)
  
  saveRDS(mgranges, file_range)
  saveRDS(mbeta, file_beta)
}

#run bsrel_matrix
bsrel_matrix('rrbs_clean_data/predictedMeth.RDS')
