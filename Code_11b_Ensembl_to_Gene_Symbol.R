# === Load Libraries ===
library(org.Hs.eg.db)

# === Define Paths ===
input_file  <- "C:\\Users\\User\\Desktop\\Soultana\\Msc Bioinformatics\\Διπλωματική-Msc\\07_Miscellaneous\\Sfikakis_Yavropoulou\\RNA-seq\\files\\Counted_Batch1_2_RawRead_PC_ranked.txt"
output_file <- "C:\\Users\\User\\Desktop\\Soultana\\Msc Bioinformatics\\Διπλωματική-Msc\\07_Miscellaneous\\Sfikakis_Yavropoulou\\RNA-seq\\files\\Counted_Batch1_2_RawRead_PC_ranked_Symbols.txt"

# === Load Ranked Matrix ===
rank_matrix <- read.table(input_file, header = TRUE, row.names = 1, sep = "\t")

# === Clean Ensembl IDs ===
# Strip version numbers from Ensembl IDs 
ensembl_ids_clean <- gsub("\\..*", "", rownames(rank_matrix))

# === Fetch Gene Symbols Locally ===
# Query the local database to map cleaned Ensembl IDs to official symbols
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys      = ensembl_ids_clean,
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "first"
)

# === Safe and Accurate Filtering ===
# 1. Create a clean mapping data frame to avoid missing values (NAs) and empty strings
mapping_df <- data.frame(
  Ensembl = rownames(rank_matrix),
  Symbol  = gene_symbols,
  stringsAsFactors = FALSE
)

# Remove rows where no official Gene Symbol was found (NAs)
mapping_df <- mapping_df[!is.na(mapping_df$Symbol), ]

# Remove rows with empty strings as Gene Symbols
mapping_df <- mapping_df[mapping_df$Symbol != "", ]

# Remove duplicate Gene Symbols to ensure unique names for row indexing
mapping_df <- mapping_df[!duplicated(mapping_df$Symbol), ]

# 2. Filter the original ranked matrix to retain only the validated genes
final_matrix <- rank_matrix[mapping_df$Ensembl, ]

# 3. Securely replace Ensembl IDs with the official unique Gene Symbols
rownames(final_matrix) <- mapping_df$Symbol

# === Save the New File ===
write.table(final_matrix, file = output_file, sep = "\t", quote = FALSE, col.names = NA)
cat("Success! The file with Gene Symbols has been created properly.\n")