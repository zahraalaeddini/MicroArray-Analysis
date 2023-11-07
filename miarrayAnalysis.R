$setwd("C:/Users/Zahra/Desktop/bio pro/")

library(GEOquery)
library(umap)
library(limma)
library(pheatmap)
library(ggplot2)
library(gplots)
library(reshape2)
library(plyr)

series <- "GSE9476"
platform <- "GPL96"

####load data
gset <- getGEO(series, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "data/")

if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

gr <- c("CD34", rep("BM", 10), rep("CD34", 7), rep("aml", 26), rep("PB", 10), rep("CD34", 10))

#### log2 scale if required
ex <- exprs(gset)
# ex <- log2(ex + 1)
# exprs(gset) <- ex

pdf("Results/boxplot.pdf", width = 64)
boxplot(ex)
dev.off()

####normalize if required
# ex <- normalizeQuantiles(ex)
# exprs(gset) <- ex

####Correlation Heatmap
pdf("Results/corHeatmap.pdf", width=15, height=15)
pheatmap:: pheatmap(cor(ex), labels_row = gr, labels_col = gr, color=bluered(256), border_color = NA)
dev.off()

####principal component analysis
pc <- prcomp(ex)
pdf("Results/PC.pdf")
plot(pc)
dev.off()

pdf("Results/PC_.pdf")
plot(pc$x[,1:2])
dev.off()

ex.scale <- t(scale(t(ex), scale=FALSE))
pc <- prcomp(ex.scale)

pdf("Results/PC_scale.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()


pcr <- data.frame(pc$r[, 1:3], Group=gr)

pdf("Results/PCA_samples.pdf")
ggplot(pcr, aes(PC1, PC2, color=Group)) + geom_point(size=3) + theme_bw()
dev.off()


#### Differential Expression Analysis
gr <- factor(gr)
gset$description <- gr
design <- model.matrix(~description + 0, gset)
colnames(design) <- levels(gr)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(aml - CD34, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by = "B", number = Inf)

tT <- subset(tT, select=c("Gene.symbol", "Gene.ID", "adj.P.Val", "logFC"))
write.table(tT, "Results/AML_Cd34.txt", row.names=F, sep="\t", quote = F)


aml.up <- subset(tT, logFC > 1 & adj.P.Val < 0.05)
#aml.up.genes <- sub("///.*", "",aml.up.genes)
aml.up.genes <- unique(as.character(strsplit2(aml.up$Gene.symbol, "///")))
write.table(aml.up.genes, file="Results/AML_CD34_UP.txt", quote = F, row.names = F, col.names = F)

aml.down <- subset(tT, logFC < -1 & adj.P.Val < 0.05)
#aml.down.genes <- sub("///.*", "", aml.down.genes)
aml.down.genes <- unique(as.character(strsplit2(aml.down$Gene.symbol, "///")))
write.table(aml.down.genes, file="Results/AML_CD34_down.txt", quote = F, row.names = F, col.names = F)



