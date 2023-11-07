setwd("C:/Users/Zahra/Desktop/batch/")

library(GEOquery)
library(umap)
library(limma)
library(pheatmap)
library(ggplot2)
library(gplots)
library(reshape2)
library(plyr)

series <- "GSE46872"
platform <- "GPL6244"
gset <- getGEO(series, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "data/")

