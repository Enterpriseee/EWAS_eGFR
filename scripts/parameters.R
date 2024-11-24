options(stringsAsFactors = FALSE)

library(openxlsx)
bm_phen <- read.xlsx("sample_information.xlsx")
samples_keep <- bm_phen$no[which(!bm_phen$no %in% c("4063", 
                                                    "3111", "3112", 
                                                    "2163", "2164", "3021", "3022", "4057", "4058", "3051", "3052", "5025", "5026",
                                                    "3005", "3006"))]
samples_keep <- as.character(samples_keep)
bm_phen <- bm_phen[bm_phen$no %in% samples_keep, ]
table(table(bm_phen$family_ID))
all(table(bm_phen$family_ID) == 2)
rm(bm_phen)

