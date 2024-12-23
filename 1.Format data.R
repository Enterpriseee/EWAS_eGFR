#library software
library(data.table)
library(openxlsx)

#create data.frame and add the contents of the sample column
b <- data.frame(
  sample = NA,
  gender = NA,
  age = NA,
  eGFR = NA
)

a<-fread("expression data 24.csv")
col_names <- colnames(a)[-1]
if (length(col_names) > nrow(b)) {
  new_rows <- length(col_names) - nrow(b)
  new_data <- data.frame(matrix(ncol = ncol(b), nrow = new_rows))
  colnames(new_data) <- colnames(b)
  b <- rbind(b, new_data)
  }
b$sample <- col_names


#add gender, age, and eGFR
library(dplyr)
a<-read.xlsx("sample_information.xlsx")
a$no <- as.character(a$no)
b$sample <- as.character(b$sample)
a$no <- paste0("S", a$no)
b <- b %>%
  left_join(a, by = c("sample" = "no")) %>%
  select(sample, gender = gender.y, age = age.y, eGFR = eGFR.y)

#add family_ID
b <- b %>%
  mutate(family_ID = NA) %>%
  select(family_ID, everything())
numbers <- rep(1:ceiling(nrow(b) / 2), each = 2)
numbers <- numbers[1:nrow(b)]
b$family_ID <- numbers

#add ENSG_id
a<-fread("ENSG.csv",header = F)
a<-unique(a)
a <- subset(a, !is.na(V1))
new_column_names <- as.character(a$V1)
a <- setNames(data.frame(matrix(ncol = length(new_column_names), nrow = 1)), new_column_names)
new_column_names <- colnames(a)
for (col_name in new_column_names) {
  b[[col_name]] <- NA
}
a<-fread("expression data 24.csv")

for (col_name in colnames(b)[-1]) {  
  for (row_index in 1:nrow(b)) {    
    gene_name <- col_name            
    sample_name <- b$sample[row_index]  
    
    gene_index <- which(a$gene == gene_name)
    if (length(gene_index) == 1) {  
      if (sample_name %in% colnames(a)) {
        b[row_index, col_name] <- a[gene_index, get(sample_name)]
      }
    }
  }
}

#save data
b <- b[!is.na(b$gender), ]
write.csv(b,file = "eGFRtraits.csv")
