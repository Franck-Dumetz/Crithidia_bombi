# This script takes a raw count table of each gene in 3 adherent and 3 swimming samples and uses DESeq2 and EdgeR
# to perform DGE analysis. This script also generates a bar graph of the number of mapped reads, a PCA plot, and
# a volcano plot.

library(DESeq2)

counts <- read.csv("counts_6-18_1.csv", header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

coldata <- data.frame(row.names = colnames(counts), condition = c("Adherent","Adherent","Adherent", "Swimming","Swimming","Swimming")
)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)

dds <- DESeq(dds)
res <- results(dds)

levels(dds$condition)

head(res)

#mapped reads plot
library(ggplot2)

df <- data.frame(
  sample = c("Adherent1", "Adherent2", "Adherent3", "Swimming1", "Swimming2", "Swimming3"),
  reads  = c(298064835, 275864051, 338993806, 213729288, 298801233, 297619740),
  group  = c("Adherent", "Adherent", "Adherent", "Swimming", "Swimming", "Swimming")
)

ggplot(df, aes(x = sample, y = reads, fill = group)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::label_scientific()) +
  scale_x_discrete(
    labels = c(
      "",
      "\nAdherent",
      "",
      "",
      "\nSwimming",
      ""
    )
  ) +
  theme_classic() +
  theme(legend.position = "none", plot.margin = margin(10, 10, 30, 10)) +
  labs(
    y = "Reads mapped to C. bombi genome",
    x = NULL,
    title = "Mapped Reads"
  )

#PCA
vsd <- vst(dds, blind = FALSE) # variance stabilization transformation (normalized counts)
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))
#Eigenvector with largest eigenvalue is PC1
ggplot(pca,aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 6) +
  theme_classic() +
  labs(
    color = NULL,
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
  )

#Volcano Plot
res$gene <- rownames(res)
res <- na.omit(res)

res$significant <- "Not significant"

res$significant[
  res$log2FoldChange >= 1 & res$padj < 0.01
] <- "Up in swimming"

res$significant[
  res$log2FoldChange <= -1 & res$padj < 0.01
] <- "Up in adherent"

library(ggplot2)

ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significant), alpha = 0.7, size = 1.5) +
  coord_cartesian(clip = "off") +
  scale_color_manual(
    values = c(
      "Up in adherent" = "red",
      "Up in swimming" = "blue",
      "Not significant" = "grey"
    )
  ) +
  theme_classic() +
  labs(
    color = NULL,
    x = "log2 fold change",
    y = "-log10 adjusted p-value",
  )


if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")


#EdgeR
library(edgeR)
counts <- read.csv("counts_6-18_1.csv", row.names=1)
group <- factor(c("Adherent","Adherent","Adherent",
                  "Swimming","Swimming","Swimming"))
dge <- DGEList(counts = counts, group = group)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)
norm <- cpm(dge)
att_mean <- rowMeans(norm[, group == "Adherent"])
swim_mean <- rowMeans(norm[, group == "Swimming"])
norm_counts_edger <- data.frame(
  Gene = rownames(norm),
  att_mean_CPM = att_mean,
  swim_mean_CPM = swim_mean
)

write.csv(norm_counts_edger, "norm_counts_bombi_6-25.csv", row.names = TRUE)
design <- model.matrix(~group)
design

dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, coef=2)
results <- topTags(lrt, n=Inf)$table
head(results)

results$gene <- rownames(results)
results <- na.omit(results)
write.csv(results, "full_edger.csv", row.names = TRUE)
sig <- results[results$PValue < 0.01 & results$logFC >= 1, ]
write.csv(sig, "up_in_swimming.csv", row.names = TRUE)

sig <- results[results$PValue < 0.01 & results$logFC <= -1, ]
write.csv(sig, "up_in_adherent.csv", row.names = TRUE)

results$significant <- "Not significant"

results$significant[
  results$logFC >= 1 & results$PValue < 0.01
] <- "Up in swimming"

results$significant[
  results$logFC <= -1 & results$PValue < 0.01
] <- "Up in adherent"

##Volcano Plot

library(ggplot2)

ggplot(results, aes(x = logFC, y = -log10(PValue))) +
  geom_point(aes(color = significant), alpha = 0.7, size = 1.5) +
  coord_cartesian(clip = "off") +
  scale_color_manual(
    values = c(
      "Up in adherent" = "red",
      "Up in swimming" = "blue",
      "Not significant" = "grey"
    )
  ) +
  theme_classic() +
  labs(
    color = NULL,
    x = "log2 fold change",
    y = "-log10 adjusted p-value",
  )

