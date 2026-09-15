# TRANSCRIPTOMIC PIPELINE: Differential Gene Expression Analysis
# Project: Cortisol-induced Stress Response 

# [1] ENVIRONMENT PREPARATION
if (!require("readxl")) install.packages("readxl")
if (!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db") # Added local annotation database
library(readxl)
library(DESeq2)
library(ggplot2)
library(org.Hs.eg.db) # Loaded for Ensembl-to-Symbol mapping

# [2] DATA PATH CONFIGURATION
# Define absolute paths for input files and output directory
excel_path  <- "C:/Users/User/Desktop/Soultana/Msc Bioinformatics/Διπλωματική-Msc/07_Miscellaneous/Sfikakis_Yavropoulou/RNA-seq/files/Low_High_Cortizol_6_4_2026.xlsx"
counts_path <- "C:/Users/User/Desktop/Soultana/Msc Bioinformatics/Διπλωματική-Msc/07_Miscellaneous/Sfikakis_Yavropoulou/RNA-seq/files/featureCounts_raw_counts.txt"
path_out    <- "C:/Users/User/Desktop/Soultana/Msc Bioinformatics/Διπλωματική-Msc/07_Miscellaneous/Sfikakis_Yavropoulou/RNA-seq/files/"

# [3] METADATA PRE-PROCESSING & FACTOR LEVEL DEFINITION
# Importing experimental design and cleaning non-biological entries
metadata <- read_excel(excel_path)
metadata <- as.data.frame(metadata)

# Filter out rows with NA values or summary statistics 
metadata <- metadata[!is.na(metadata[[1]]), ] 
metadata <- metadata[metadata[[1]] != "AVG" & metadata[[1]] != "TTEST", ]
colnames(metadata)[1] <- "sample_id"

# Fix duplicates
metadata <- metadata[!duplicated(metadata$sample_id), ]

metadata$condition <- NA
metadata$condition[1:20] <- "high"
metadata$condition[21:nrow(metadata)] <- "low"

# [4] TRANSCRIPT COUNT MATRIX IMPORT & CLEANING
# Load raw read counts generated from featureCounts
counts <- read.table(counts_path, header = TRUE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)

# Fix duplicates in IDs
counts <- counts[!duplicated(counts[,1]), ]

# Set Ensembl IDs as row identifiers
rownames(counts) <- counts[,1]
counts_matrix <- counts[, -(1)]

# Cast to integer matrix: Essential requirement for Negative Binomial distribution models
counts_matrix <- as.matrix(counts_matrix)
mode(counts_matrix) <- "integer"

# Regex transformation to standardize sample names
sample_codes <- sub(".*(25K[0-9]{3}).*", "\\1", colnames(counts_matrix))
colnames(counts_matrix) <- sample_codes

# Replicate aggregation 
if(any(duplicated(colnames(counts_matrix)))) { 
message("Duplicates found in counts columns. Aggregating...") 
counts_matrix <- t(rowsum(t(counts_matrix), group = colnames(counts_matrix))) }

# [5] BIOINFORMATIC DATA ALIGNMENT & INTEGRITY CHECKS
# Identify intersection of samples present in both counts and metadata
common_samples <- intersect(colnames(counts_matrix), metadata$sample_id)

# Subset both objects to ensure 1:1 mapping
counts_matrix <- counts_matrix[, common_samples]
metadata      <- metadata[metadata$sample_id %in% common_samples, ]

# Remove any rows where condition is still NA 
metadata <- metadata[!is.na(metadata$condition), ]

# Reorder metadata rows to match count matrix column order 
metadata <- metadata[match(colnames(counts_matrix), metadata$sample_id), ] 
rownames(metadata) <- metadata$sample_id

metadata$condition <- factor(metadata$condition, levels = c("low", "high"))

# Export synchronized raw count matrix for downstream reproducibility
write.table(counts_matrix, 
            file = paste0(path_out, "featureCounts_raw_counts_final.bed"), 
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# [6] STATISTICAL ANALYSIS: DESeq2 Pipeline
# Constructing DESeqDataSet object from matrix and colData
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = metadata,
  design = ~ condition
)

# Estimate size factors, dispersions, and execute Wald test
dds <- DESeq(dds)

# [7] RESULTS EXTRACTION & SIGNIFICANCE FILTERING
# Contrasting High vs Low Cortisol expression profiles
res <- results(dds, contrast = c("condition", "high", "low"))
res_df <- as.data.frame(res)

# --- Start of Ensembl to Gene Symbol Mapping ---
# Clean Ensembl IDs by removing version numbers 
ensembl_clean <- gsub("\\..*", "", rownames(res_df))

# Fetch symbols locally
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys      = ensembl_clean,
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "first"
)

# Build a safe mapping framework
mapping_df <- data.frame(
  Ensembl = rownames(res_df),
  Symbol  = gene_symbols,
  stringsAsFactors = FALSE
)

# Filter out rows with missing or empty symbols, and drop duplicates
mapping_df <- mapping_df[!is.na(mapping_df$Symbol), ]
mapping_df <- mapping_df[mapping_df$Symbol != "", ]
mapping_df <- mapping_df[!duplicated(mapping_df$Symbol), ]

# Subset and update the statistical dataframe to use official Gene Symbols
res_df <- res_df[mapping_df$Ensembl, ]
rownames(res_df) <- mapping_df$Symbol
# --- End of Mapping ---

# Extract Significant Differentially Expressed Genes 
sig_genes <- res_df[which(res_df$padj < 0.05), ]

# Save comprehensive results and significant subset to TSV format
write.table(res_df, 
            file = paste0(path_out, "DESeq_results.tsv"), 
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

write.table(sig_genes, 
            file = paste0(path_out, "DESeq_results_sig.tsv"), 
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)



# DIFFERENTIAL EXPRESSION CATEGORIZATION 
# Initialize a classification vector for significant and non-significant features
res_df$diffexpressed <- "Non-Significant"

# Define Up-regulated genes (Overexpressed): 
# Criteria: Log2 Fold Change >= 1 (at least 2-fold increase) AND p-value <= 0.05
res_df$diffexpressed[res_df$log2FoldChange >= 1 & res_df$pvalue <= 0.05] <- "Overexpressed"

# Define Down-regulated genes (Underexpressed): 
# Criteria: Log2 Fold Change <= -1 (at least 2-fold decrease) AND p-value <= 0.05
res_df$diffexpressed[res_df$log2FoldChange <= -1 & res_df$pvalue <= 0.05] <- "Underexpressed"

# Refactor labels into a categorical variable for consistent plotting legend ordering
res_df$diffexpressed <- factor(res_df$diffexpressed, levels = c("Overexpressed", "Underexpressed", "Non-Significant"))


# VOLCANO PLOT VISUALIZATION 
if(!require(ggrepel)) install.packages("ggrepel")
library(ggrepel)
res_plot <- as.data.frame(res_df)
res_plot <- res_plot[!is.na(res_plot$pvalue), ]

# Mapping genomic features: X-axis represents magnitude of change, Y-axis represents statistical confidence
ggplot(data = res_df, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed)) +
  
# Scatter plot implementation with alpha blending to mitigate overplotting
geom_point(alpha = 0.4, size = 1.5) + 
theme_minimal() +
  
# Custom aesthetic scale for distinct biological status representation
scale_color_manual(values = c("Overexpressed" = "red", 
                              "Underexpressed" = "blue", 
                              "Non-Significant" = "grey")) +
  
geom_text_repel(data = subset(res_plot, diffexpressed != "Non-Significant"), 
                 aes(label = rownames(subset(res_plot, diffexpressed != "Non-Significant"))),
                 size = 2.5, max.overlaps = 10, show.legend = FALSE) +
  
# Add biological and statistical significance thresholds (Horizontal: p=0.05, Vertical: 2nd fold)
geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed") +
geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed") +
theme_minimal() +
  
# Annotate plot axes and descriptive titles
labs(title = "Volcano Plot: Differential Gene Expression Profile",
     x = "Log2 Fold Change",
     y = "-Log10 P-value",
     color = "Biological Status")




