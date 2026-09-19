#!/usr/bin/env python3
# === MATRIX for Top 500 Genes ===
# Project: Cortisol-induced Stress Response Transcriptomic Analysis
# Description: This script filters the full transcriptomic ranked matrix to retain 
#              only the top 500 most important genes identified by the Random Forest classifier.

import pandas as pd

# === Define Absolute File Paths ===
# 1. Input: The original full matrix with all 78,714 genes across the 40 biological samples
full_matrix_path = r"C:\Users\User\Desktop\Soultana\Msc Bioinformatics\Διπλωματική-Msc\07_Miscellaneous\Sfikakis_Yavropoulou\RNA-seq\files\Counted_Batch1_2_RawRead_PC_ranked.txt"

# 2. Input: The finalized table of the top 500 genes (with Ensembl IDs and Gene Symbols)
top_500_mapped_path = r"C:\Users\User\Desktop\Soultana\Msc Bioinformatics\Διπλωματική-Msc\07_Miscellaneous\Sfikakis_Yavropoulou\RNA-seq\files\Final_Table_Top_500_Genes.txt"

# 3. Output: The target smaller expression matrix for the top 500 features
output_matrix_path = r"C:\Users\User\Desktop\Soultana\Msc Bioinformatics\Διπλωματική-Msc\07_Miscellaneous\Sfikakis_Yavropoulou\RNA-seq\files\Filtered_Top500_Ranked_Matrix.txt"

# Load the full genomic matrix 
full_df = pd.read_csv(full_matrix_path, sep='\t')

# Ensure the first column is named 'Ensembl_ID' to perform a flawless cross-reference mapping
full_df.rename(columns={full_df.columns[0]: 'Ensembl_ID'}, inplace=True)

# Load the top 500 annotated features table
top_500_df = pd.read_csv(top_500_mapped_path, sep='\t')

# This operation automatically discards the remaining 78,214 non-top features
merged_df = pd.merge(top_500_df[['Ensembl_ID']], full_df, on='Ensembl_ID', how='inner')

# === Export Finalized Subset Matrix ===
merged_df.to_csv(output_matrix_path, sep='\t', index=False)

print(f"\nSuccess! The target subset matrix was successfully generated.")
print(f"File saved to: {output_matrix_path}")
print(f"Matrix Dimensions -> Rows (Genes): {merged_df.shape[0]} | Columns (Ensembl ID + 40 Samples): {merged_df.shape[1]}")





import pandas as pd

# 1. Ορισμός διαδρομής και φόρτωση καθολικού πίνακα
full_matrix_path = r"C:\Users\User\Desktop\Soultana\Msc Bioinformatics\Διπλωματική-Msc\07_Miscellaneous\Sfikakis_Yavropoulou\RNA-seq\files\Counted_Batch1_2_RawRead_PC_ranked.txt"
full_df = pd.read_csv(full_matrix_path, sep="\t")
full_df.rename(columns={full_df.columns[0]: "Ensembl_ID"}, inplace=True)

# 2. Δημιουργία του X_full (μετά το transpose τα δείγματα είναι στις γραμμές)
X_full = full_df.drop(columns=["Ensembl_ID"]).T

# 3. Εκτύπωση των ονομάτων των δειγμάτων στη σειρά
print("Σειρά δειγμάτων στον πίνακα:")
for i, sample_name in enumerate(X_full.index):
  print(f"{i}: {sample_name}")