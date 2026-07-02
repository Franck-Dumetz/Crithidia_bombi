# This script uses previous a list of upregulated C.bombi genes, finds the corresponding orthologs and gene
# ontology terms, and performed GO enrichment analysis using topGO. It generates plots showing how many genes
# correspond to each of the top GO terms and fold enrichment of significant genes.

library(edgeR)
counts <- read.csv("Downloads/counts-615.csv", row.names=1)
group <- factor(c("Adherent","Adherent","Adherent",
                  "Swimming","Swimming","Swimming"))
dge <- DGEList(counts = counts, group = group)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)
design <- model.matrix(~group)
design

dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, coef=2)
results <- topTags(lrt, n=Inf)$table

install.packages("BiocManager", repos = "http://cran.us.r-project.org")
BiocManager::install("topGO")
library(topGO)

up <- read.csv("up_in_adherent.csv", row.names=1)
up_orthologs <- ifelse(grepl("cf_ortholog=", up$gene), sub(".*cf_ortholog=([^;]+).*", "\\1", up$gene), NA)
keep <- !is.na(up_orthologs) & up_orthologs != ""
up_orthologs <- up_orthologs[keep]
up_orthologs <- sub("\\.1$", "", up_orthologs)
up_orthologs <- unique(up_orthologs)

edger <- read.csv("full_edger.csv", row.names=1)
edger_orthologs <- ifelse(grepl("cf_ortholog=", edger$gene), sub(".*cf_ortholog=([^;]+).*", "\\1", edger$gene), NA)
keep <- !is.na(edger_orthologs) & edger_orthologs != ""
edger_orthologs <- edger_orthologs[keep]
edger_orthologs <- sub("\\.1$", "", edger_orthologs)
edger_orthologs <- unique(edger_orthologs)

gene_list <- factor(as.integer(edger_orthologs %in% up_orthologs))
names(gene_list) <- edger_orthologs

gaf <- read.delim("~/Downloads/TriTrypDB-68_CfasciculataCfCl_Curated_GO.gaf", comment.char = "!", header = FALSE, stringsAsFactors = FALSE)
gene2go <- split(gaf$V5, gaf$V2)
sum(names(gene_list) %in% names(gene2go))

GOdata <- new("topGOdata", ontology = "MF", allGenes = gene_list, annot = annFUN.gene2GO, gene2GO = gene2go)
results <- runTest(GOdata, algorithm = "weight01", statistic = "Fisher")
allGO <- usedGO(GOdata)
table <- GenTable(GOdata, results, orderBy="results", topNodes=20)
table

table$Significant <- as.numeric(table$Significant)
table$Expected <- as.numeric(table$Expected)

table$FoldEnrichment <- table$Significant / table$Expected
table <- table[as.numeric(table$result1) < 0.05, ]
top <- table$GO.ID
gene_list_by_GO <- lapply(top, function(go) {
  genes <- genesInTerm(GOdata, go)[[1]]
  genes[genes %in% names(gene_list)[gene_list == 1]]
})
names(gene_list_by_GO) <- top
gene_list_by_GO

table

#write.table(table, "Downloads/GO-files/GO_up_swim_BP.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#capture.output(gene_list_by_GO, file = "Downloads/GO-files/GO_up_swim_BP_genes.txt")

# Get number of C.bombi genes

cbombi_to_ortholog <- data.frame(
  cbombi_gene = rownames(up),
  ortholog = ifelse(grepl("cf_ortholog=", up$gene),
                    sub(".*cf_ortholog=([^;]+).*", "\\1", up$gene),
                    NA)
)

cbombi_to_ortholog$ortholog <- sub("\\.1$", "", cbombi_to_ortholog$ortholog)
cbombi_to_ortholog <- cbombi_to_ortholog[!is.na(cbombi_to_ortholog$ortholog), ]
cbombi_counts <- sapply(top, function(go) {
  orthologs <- genesInTerm(GOdata, go)[[1]]
  
  sum(cbombi_to_ortholog$ortholog %in% orthologs)
})

library(ggplot2)
table <- GenTable(GOdata, results, orderBy="results", topNodes=20)
table$Significant <- as.numeric(table$Significant)
cbombi_counts <- sapply(table$GO.ID, function(go) {
  orthologs <- genesInTerm(GOdata, go)[[1]]
  sum(cbombi_to_ortholog$ortholog %in% orthologs)
})
table$Count <- cbombi_counts

ggplot(table,
       aes(x = Count,
           y = reorder(Term, Count),
           fill = Count)) +
  geom_col() +
  labs(
    x = "Number of Genes",
    y = NULL,
    fill = "Gene Count"
  )  +
  scale_fill_gradient(
    low = "black",
    high = "#d95f02"
  ) +
  theme_bw()

ggplot(table,
       aes(x = FoldEnrichment,
           y = reorder(Term, FoldEnrichment),
           fill = FoldEnrichment)) +
  geom_col() +
  theme_bw() +
  labs(
    x = "Fold Enrichment (Observed / Expected)",
    y = NULL,
    fill = "Fold Enrichment"
  )

# Enrichment plot
table$result1 <- as.numeric(table$result1)
table$score <- -log10(table$result1)

ggplot(table,
       aes(x = reorder(Term, score),
           y = score)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "GO Biological Process Enrichment",
    x = "GO Term",
    y = "-log10(p-value)"
  ) +
  theme_bw()

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GSEABase", "GO.db", "AnnotationDbi"))

library(GSEABase)
library(GO.db)
library(AnnotationDbi)
#genes <- readLines("gene_ids_trimmed.txt")
#all_genes <- readLines("all_gene_ids.txt")

#gene_list <- factor(as.integer(all_genes %in% genes))
#names(gene_list) <- all_genes

#mappings <- readMappings(file = "mapping_ids.txt")

#GOdata <- new("topGOdata", ontology = "MF", allGenes = gene_list, annot = annFUN.gene2GO, gene2GO = mappings)

#results <- runTest(GOdata, algorithm = "weight01", statistic = "Fisher")

#table <- GenTable(GOdata, results, topNodes = 20)

#write.table(table, "GO_enrichment_transcripts.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#sig_genes <- genesInTerm(GOdata, "GO:0008270")[[1]]
#sig_genes <- intersect(sig_genes, genes)

#print(sig_genes)



