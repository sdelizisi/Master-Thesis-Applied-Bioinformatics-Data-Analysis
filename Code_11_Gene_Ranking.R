#!/usr/bin/env Rscript

# Code 11
# === Transform Counts to Ranks ===

# === Define Input and Output Files ===
input_file  <- "C:\\Users\\User\\Desktop\\Soultana\\Msc Bioinformatics\\Διπλωματική-Msc\\07_Miscellaneous\\Sfikakis_Yavropoulou\\RNA-seq\\files\\GeTMM_normalized_counts.txt"
output_file <- "C:\\Users\\User\\Desktop\\Soultana\\Msc Bioinformatics\\Διπλωματική-Msc\\07_Miscellaneous\\Sfikakis_Yavropoulou\\RNA-seq\\files\\Counted_Batch1_2_RawRead_PC_ranked.txt"

# === Load Input Matrix ===
# Load the input matrix from the specified file (genes as rows, samples as columns)
matrix <- read.table(input_file, header = TRUE, row.names = 1, sep = "\t")
cat("Input matrix loaded successfully.\n")

# === Rank Genes Within Each Sample ===
# Apply the 'rank' function to each column (sample), resolving ties by assigning the minimum rank
rank_matrix <- apply(matrix, 2, rank, ties.method = "min")
cat("Gene ranking completed.\n")

# === Write Output Matrix ===
# Write the ranked matrix to the specified output file
write.table(rank_matrix, file = output_file, sep = "\t", quote = FALSE, col.names = NA)
cat("Ranked matrix written to file successfully!\n")
