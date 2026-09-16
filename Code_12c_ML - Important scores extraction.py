#!/usr/bin/env python3
# === ML Feature Importance Extraction (All and Top 500) ===
# Project: Cortisol-induced Stress Response Transcriptomic Analysis

import random
from numpy import genfromtxt
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier

# === Seed Initialization for Reproducibility ===
np.random.seed(1)
random.seed("Soultana")


# === Load Data Function ===
def load_data(dataset_path):
  dataset = genfromtxt(dataset_path, delimiter="\t", dtype="str")

  # Extract Ensembl Gene IDs from the first column (excluding header)
  gene_ids = dataset[1:, 0]

  # Extract numeric matrix and transpose (samples as rows, genes as columns)
  data = dataset[1:, 1:].astype("float").T

  print(f"Data Shape verified (Samples, Genes): {data.shape}")
  print(f"Total Gene IDs loaded: {len(gene_ids)}")
  return data, gene_ids


# === Load Labels Function ===
def load_labels():
  # Balanced cohort: 20 High (1) and 20 Low (0) Cortisol samples
  labels = [1] * 20 + [0] * 20
  print(f"Target labels loaded. Total sample count: {len(labels)}")
  return np.array(labels)


# === Feature Importance Extraction Function ===
def extract_feature_importance(
    X, y, gene_ids, output_all_path, output_top500_path
):
  classifier = RandomForestClassifier(criterion="gini", random_state=1)
  classifier.fit(X, y)

  importances = classifier.feature_importances_

  # Alignment verification between features and gene IDs
  if len(gene_ids) != len(importances):
    print(
        f"Alignment warning: Gene IDs ({len(gene_ids)}) and model features"
        f" ({len(importances)}) mismatch. Aligning..."
    )
    gene_ids = gene_ids[: len(importances)]

  # Create DataFrame and sort in descending order of importance
  importance_df = pd.DataFrame(
      {"Ensembl_ID": gene_ids, "Importance_Score": importances}
  ).sort_values(by="Importance_Score", ascending=False)

  # 1. Export the full table with ALL genes (78,724)
  importance_df.to_csv(output_all_path, sep="\t", index=False)
  print(f"1. Complete dataset (All genes) saved to:\n   {output_all_path}")

  # 2. Isolate and export the top 500 most influential gene features
  top_500_genes = importance_df.head(500)
  top_500_genes.to_csv(output_top500_path, sep="\t", index=False)
  print(f"2. Top 500 features saved to:\n   {output_top500_path}")


# === Main Function ===
def main():
  ranked_matrix_path = r"C:\Users\User\Desktop\Soultana\Msc Bioinformatics\Διπλωματική-Msc\07_Miscellaneous\Sfikakis_Yavropoulou\RNA-seq\files\Counted_Batch1_2_RawRead_PC_ranked.txt"
  output_all_path = r"C:\Users\User\Desktop\Soultana\Msc Bioinformatics\Διπλωματική-Msc\07_Miscellaneous\Sfikakis_Yavropoulou\RNA-seq\files\all_importance_genes.txt"
  output_top500_path = r"C:\Users\User\Desktop\Soultana\Msc Bioinformatics\Διπλωματική-Msc\07_Miscellaneous\Sfikakis_Yavropoulou\RNA-seq\files\top_500_importance_genes.txt"

  # Load data and labels (single pass read)
  data, gene_ids = load_data(ranked_matrix_path)
  labels = load_labels()

  # Train on 100% of data (N=40) and extract importance rankings
  extract_feature_importance(
      data, labels, gene_ids, output_all_path, output_top500_path
  )


if __name__ == "__main__":
  main()