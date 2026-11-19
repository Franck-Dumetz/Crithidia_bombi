# Differential protein abundance, C. bombi WHA1 attached vs swimming
DATA    <- Sys.getenv("CB_DATA",    "../data")
RESULTS <- Sys.getenv("CB_RESULTS", "../results")

# limma on WHA1 MaxLFQ, excluding the 566 Cb.1.* entries (08076 gene set)
suppressPackageStartupMessages({library(limma)})
d <- read.delim(file.path(DATA, "prot_raw.tsv"), check.names=FALSE)
d <- d[!grepl("__CONT", d$Description), ]
d$pid <- trimws(as.character(d$`Protein ID`))
n0 <- nrow(d)
d <- d[grepl("^CbWHA1", d$pid), ]
cat(sprintf("rows: %d non-contaminant -> %d after dropping Cb.1.* (%d removed)\n",
            n0, nrow(d), n0 - nrow(d)))
S <- c(paste0("WHA1Attached-",1:3), paste0("WHA1Swimming-",1:3))
M <- as.matrix(d[, paste(S,"MaxLFQ Intensity")]); colnames(M) <- S
M[M==0] <- NA; L <- log2(M); rownames(L) <- trimws(as.character(d$Protein))
A <- S[1:3]; W <- S[4:6]
nA <- rowSums(!is.na(L[,A])); nW <- rowSums(!is.na(L[,W]))
keep <- d$`Combined Total Peptides` >= 3 & nA >= 2 & nW >= 2
cat(sprintf("protein groups tested: %d\n", sum(keep)))
E <- L[keep,]
grp <- factor(c(rep("att",3), rep("swim",3)), levels=c("swim","att"))
fit <- eBayes(lmFit(E, model.matrix(~grp)), trend=TRUE, robust=TRUE)
tt <- topTable(fit, coef=2, number=Inf, sort.by="none")
g <- sub("\\.[0-9]+(\\.p[0-9]+)?:.*$", "", rownames(E), perl=TRUE)
  out <- data.frame(orf=rownames(E), protein=g, prot_logFC=tt$logFC, prot_P=tt$P.Value,
                  prot_FDR=tt$adj.P.Val, prot_AveExpr=tt$AveExpr)
write.table(out, file.path(RESULTS, "prot_limma.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf("significant (FDR<0.05 & |logFC|>1): %d  (up att %d, up swim %d)\n",
    sum(tt$adj.P.Val<0.05 & abs(tt$logFC)>1),
    sum(tt$adj.P.Val<0.05 & tt$logFC>1), sum(tt$adj.P.Val<0.05 & tt$logFC< -1)))
