# This script takes a raw count table of each gene in 3 adherent and 3 swimming samples and uses DESeq2 and EdgeR
# to perform DGE analysis. This script also generates a PCA plot, and a volcano plot.

library(edgeR)
counts <- read.csv("norm_counts_7-28.csv", row.names=1)
cols <- c("att1", "att2", "att3", "swim1", "swim2", "swim3")
counts_filt <- counts[rowSums(counts[cols] > 10) >= 3, ]
counts <- counts_filt

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

write.csv(norm_counts_edger, "norm_counts_bombi_7-30.csv", row.names = TRUE)

norm_counts <- read.csv("norm_counts_7-28.csv", row.names = 1)
group <- factor(c("Adherent","Adherent","Adherent",
                  "Swimming","Swimming","Swimming"))
#att_mean <- rowMeans(norm_counts[, group == "Adherent"])
#swim_mean <- rowMeans(norm_counts[, group == "Swimming"])

identical(names(swim_mean), names(att_mean))
fc = swim_mean/att_mean

#att_max = apply(norm_counts[, group == "Adherent"], 1, max)
#swim_max = apply(norm_counts[, group == "Swimming"], 1, max)

design <- model.matrix(~group)
design

dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, coef=2)
results <- topTags(lrt, n=Inf)$table
head(results)

results$gene <- rownames(results)
results <- na.omit(results)

fc <- fc[rownames(results)]
#swim_max <- swim_max[rownames(results)]
#att_max <- att_max[rownames(results)]

#identical(rownames(results), names(att_max))
results$FC <- log2(fc)
#results$AttMax <- att_max
#results$SwimMax <- swim_max
write.csv(results, "full_edger.csv", row.names = TRUE)
sig <- results[results$FDR < 0.01 & results$logFC >= 1,]
write.csv(sig, "up_in_swimming.csv", row.names = TRUE)

sig <- results[results$FDR < 0.01 & results$logFC <= -1, ]
write.csv(sig, "up_in_adherent.csv", row.names = TRUE)

results$significant <- "Not significant"

results$significant[
  results$logFC >= 1 & results$FDR < 0.01
] <- "Up in swimming"

results$significant[
  results$logFC <= -1 & results$FDR < 0.01
] <- "Up in adherent"

##Volcano Plot

library(ggplot2)

ggplot(results, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = significant), alpha = 0.7, size = 1.5) +
  coord_cartesian(clip = "off") +
  scale_color_manual(
    values = c(
      "Up in adherent" = "#D55E00",
      "Up in swimming" = "#0072B2",
      "Not significant" = "grey"
    )
  ) +
  theme_classic() +
  labs(
    color = NULL,
    x = "log2 fold change",
    y = "-log10 FDR",
  )

logCPM <- cpm(dge, log = TRUE, prior.count = 1)
pca <- prcomp(t(logCPM), scale = FALSE)

# Percent variance explained
percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)))

# Make dataframe for ggplot
pcaData <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  condition = group,
  name = colnames(logCPM)
)

ggplot(pcaData, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 6) +
  scale_color_manual(values = c("Adherent" = "#D55E00", 
                                "Swimming" = "#0072B2")) +
  theme_classic() +
  labs(
    color = NULL,
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance")
  )
