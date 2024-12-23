#library software
options(stringsAsFactors = FALSE)
library(parallel)
library(BiSeq)
library(openxlsx)

#load data
bm_files <- list.files("methylation_data", full.names = TRUE)
bm_ids <- sub(".*([0-9]{4}).*", "\\1", bm_files, perl = TRUE)
bm_phen <- read.xlsx("sample_information.xlsx")
rrbs <- readBismark(bm_files, colData = DataFrame(row.names = bm_ids))

#Create an output directory and save the data
dir.create("rrbs_clean_data")
saveRDS(rrbs, "rrbs_clean_data/rrbs.RDS")

#Plot the coverage box plot of methylation data
png("rrbs_clean_data/rrbs_covBoxplots.png", width = 900, height = 480, units = 'px', pointsize = 12)
covBoxplots(rrbs, col = "cornflowerblue", las = 2)
dev.off()

#Cluster analysis
rrbs.clust.unlim <- clusterSites(object = rrbs,
                                 perc.samples = 4/5,
                                 min.sites = 20,
                                 max.dist = 100,
                                 mc.cores = 1)

#Filter low coverage data
ind.cov <- totalReads(rrbs.clust.unlim) > 0
quant <- quantile(totalReads(rrbs.clust.unlim)[ind.cov], 0.9)
quant
rrbs.clust.lim <- limitCov(rrbs.clust.unlim, maxCov = quant)
saveRDS(rrbs.clust.lim, 'rrbs_clean_data/rrbs.clust.lim.RDS')

#Plot the coverage box plot of the filtered data
png('rrbs_clean_data/rrbs.clust.lim_covBoxplots.png', width = 900, height = 480, units = 'px', pointsize = 12)
covBoxplots(rrbs.clust.lim, col = "cornflowerblue", las = 2)
dev.off()

#Predict methylation levels and save data
predictedMeth <- predictMeth(object = rrbs.clust.lim, mc.cores = 1)
saveRDS(predictedMeth, 'rrbs_clean_data/predictedMeth.RDS')
