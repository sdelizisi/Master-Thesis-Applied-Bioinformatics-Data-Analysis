#Code 10

# === GeTMM Normalization Script ===

import pandas as pd
import numpy as np

# 1. Load Data
count_matrix_file = r"C:\Users\User\Desktop\Soultana\Msc Bioinformatics\Διπλωματική-Msc\07_Miscellaneous\Sfikakis_Yavropoulou\RNA-seq\files\featureCounts_raw_counts.txt"
annotation_file = r"C:\Users\User\Desktop\Soultana\Msc Bioinformatics\Διπλωματική-Msc\07_Miscellaneous\Sfikakis_Yavropoulou\RNA-seq\files\featureCounts_annotation.txt"

count_matrix_df = pd.read_csv(count_matrix_file, sep="\t", index_col=0)
gene_df = pd.read_csv(annotation_file, sep="\t", index_col=0)

# 2. Gene Alignment & Length Extraction (in Kilobases)
common_genes = count_matrix_df.index.intersection(gene_df.index)
count_matrix_df = count_matrix_df.loc[common_genes]
gene_lengths_kb = gene_df.loc[common_genes, "Length"] / 1000.0

# 3. Step 1 of GeTMM: Divide by Gene Length (RPK)
rpk_df = count_matrix_df.div(gene_lengths_kb, axis=0)

# 4. Step 2 of GeTMM: Standard TMM Scale on RPK matrix
rpk_transposed = rpk_df.transpose()
rpk_lib_sizes = rpk_transposed.sum(axis=1)

# Safe log-proportions method (prevents inf/underflow)
log_rpk = np.log(rpk_transposed.div(rpk_lib_sizes, axis=0) + 1e-8)
reference_sample = log_rpk.mean(axis=0)
m_values = log_rpk.sub(reference_sample, axis=1)

# Trim top and bottom 30% of M-values to remove extreme outliers
trimmed_m = m_values.apply(lambda row: row[(row >= row.quantile(0.3)) & (row <= row.quantile(0.7))].mean(), axis=1).fillna(0.0)
geTMM_factors = np.exp(trimmed_m)

# 5. Final Scaling & Log Transformation
rpk_scaled = rpk_transposed.div(geTMM_factors, axis=0)
getmm_matrix = rpk_scaled.div(rpk_scaled.sum(axis=1), axis=0) * 1e6

# Apply log2(GeTMM + 1) for Variance Stabilization
final_getmm_log2 = np.log2(getmm_matrix + 1)
final_getmm_output = final_getmm_log2.transpose()

# 6. Save Output
output_geTMM_file = r"C:\Users\User\Desktop\Soultana\Msc Bioinformatics\Διπλωματική-Msc\07_Miscellaneous\Sfikakis_Yavropoulou\RNA-seq\files\GeTMM_normalized_counts.txt"
final_getmm_output.to_csv(output_geTMM_file, sep="\t")
print("The file GeTMM_normalized_counts.txt was successfully created!")


