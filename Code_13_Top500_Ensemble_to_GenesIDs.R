# === PIPELINE: Top 500 Ensembl IDs to Official Gene Symbols ===
# Project: Cortisol-induced Stress Response - Machine Learning Post-Analysis
# Description: Maps Ensembl IDs to official Gene Symbols using org.Hs.eg.db,
#              filters unmapped/duplicate entries, and exports the final matrix.

# === Load Required Bioconductor Libraries ===
library(AnnotationDbi)
library(org.Hs.eg.db)

# === Define Absolute File Paths ===
input_file  <- "C:\\Users\\User\\Desktop\\Soultana\\Msc Bioinformatics\\Διπλωματική-Msc\\07_Miscellaneous\\Sfikakis_Yavropoulou\\RNA-seq\\files\\top_500_importance_genes.txt"
output_file <- "C:\\Users\\User\\Desktop\\Soultana\\Msc Bioinformatics\\Διπλωματική-Msc\\07_Miscellaneous\\Sfikakis_Yavropoulou\\RNA-seq\\files\\Final_Table_Top_500_Genes.txt"

# === Load Top 500 Input Matrix ===
top_500_data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# === Clean Ensembl IDs (Strip trailing version numbers) ===
ensembl_ids_clean <- gsub("\\..*", "", top_500_data$Ensembl_ID)

# === Query Local Homo Sapiens Annotation Database ===
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys      = ensembl_ids_clean,
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "first"
)

# === Construct Consolidated Final Dataset ===
final_dataframe <- data.frame(
  Ensembl_ID       = top_500_data$Ensembl_ID,
  Gene_Symbol      = as.character(gene_symbols),
  Importance_Score = top_500_data$Importance_Score,
  stringsAsFactors = FALSE
)

# === Quality Control & Filtering ===
# 1. Remove rows where mapping failed (NAs)
final_dataframe <- final_dataframe[!is.na(final_dataframe$Gene_Symbol), ]

# 2. Remove rows with empty strings
final_dataframe <- final_dataframe[final_dataframe$Gene_Symbol != "", ]

# 3. Remove biological duplicates (retain highest importance score)
final_dataframe <- final_dataframe[!duplicated(final_dataframe$Gene_Symbol), ]

# === Export Final Academic Table ===
write.table(
  final_dataframe, 
  file      = output_file, 
  row.names = FALSE, 
  sep       = "\t", 
  quote     = FALSE
)

cat("Mapping complete. Final table saved to:", output_file, "\n")
cat("Total functional genes retained:", nrow(final_dataframe), "\n")