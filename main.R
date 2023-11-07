setwd("C:/Users/Zahra/Desktop/bio pro/data/")

library(data.table)
library(limma)

data <- read.delim("GSE9476_series_matrix.txt", comment.char = "!")
data[, -1] <- log2(1+data[, -1])
annot <- fread("GPL96.soft", skip = "!platform_table_begin", data.table = F)

colnames(data) <- c("NAME", paste0("CD34_", 1), paste0("BM_", 1:10), paste0("CD34_", 2:8), paste0("aml_", 1:26), paste0("PB_", 1:10), paste0("CD34_", 9:18))
write.table(data, "C:/Users/Zahra/Desktop/AML_GSEA.txt", row.names = F, col.names = T, sep = "\t", quote = F)
