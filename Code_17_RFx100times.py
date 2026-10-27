#!/usr/bin/env python3
# === Random Forest Robustness & AUC Evaluation (100 Iterations) ===
# This script runs Random Forest 100 times, performing a unique 90/10 split 
# at each iteration. It saves the 100 feature importance files AND computes 
# the AUC score for each run to evaluate both feature and prediction stability.

import os
import sys
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc

# === Load Data Function ===
def load_data(dataset_path):
    print(f"Loading ranked expression matrix from: {dataset_path}")
    df = pd.read_csv(dataset_path, sep='\t', index_col=0)
    feature_names = df.index.astype(str).str.strip().tolist()
    data = df.values.T
    print(f"--> Data Layout Loaded: {data.shape}")
    return data, feature_names

# === Load Labels Function ===
def load_labels():
    labels = [1] * 20 + [0] * 20
    return np.array(labels)

# === Main Pipeline ===
def main():
    # Define primary directory path
    base_dir = r"C:\Users\User\Desktop\Soultana\Msc Bioinformatics\Διπλωματική-Msc\07_Miscellaneous\Sfikakis_Yavropoulou\RNA-seq\files"
    ranked_matrix_path = os.path.join(base_dir, "Counted_Batch1_2_RawRead_PC_ranked_Symbols.txt")
    
    # Create a subfolder for the 100 files to keep the main folder clean
    output_subfolder = os.path.join(base_dir, "RF_100_Runs_Robustness_with_AUC")
    if not os.path.exists(output_subfolder):
        os.makedirs(output_subfolder)
        print(f"Created dedicated subfolder: {output_subfolder}")
    
    # 1. Load data
    data, feature_names = load_data(ranked_matrix_path)
    labels = load_labels()
    
    # List to store AUC scores from all 100 runs
    all_auc_scores = []
    
    print("\n Generating 100 individual files and calculating AUC scores. Please wait...")

    # 2. Iterative execution loop (1 to 100)
    for iteration in range(1, 101):
        # Perform a dynamic 90/10 split changing the random_state at each run
        X_train, X_test, y_train, y_test = train_test_split(
            data, labels, test_size=0.1, random_state=iteration, stratify=labels
        )
        
        # Fit the classifier using the current iteration count as the random seed
        clf = RandomForestClassifier(criterion='gini', random_state=iteration, n_jobs=-1)
        clf.fit(X_train, y_train)
        
        # Calculate AUC for this specific run on the test set
        probs = clf.predict_proba(X_test)[:, 1]
        fpr, tpr, _ = roc_curve(y_test, probs)
        run_auc = auc(fpr, tpr)
        all_auc_scores.append(run_auc)
        
        # Build a temporary dataframe for the feature importances
        df_current_run = pd.DataFrame({
            'Gene_Symbol': feature_names,
            'Importance_Score': clf.feature_importances_
        })
        
        # Sort features by importance descending (highest weights first)
        df_current_run = df_current_run.sort_values(by='Importance_Score', ascending=False)
        
        # Construct specific filename containing both run ID and its AUC score
        filename = f"RF_run_{iteration}_AUC_{run_auc:.4f}_importance.txt"
        file_export_path = os.path.join(output_subfolder, filename)
        
        # Save current run as tab-separated file
        df_current_run.to_csv(file_export_path, sep='\t', index=False)
        

    # 3. Print Final Summary Metrics
    mean_auc = np.mean(all_auc_scores)
    std_auc = np.std(all_auc_scores)
    
    print("\n" + "="*50)
    print(">>> FINAL ROBUSTNESS & PERFORMANCE METRICS <<<")
    print(f"Mean AUC Score across 100 runs: {mean_auc:.4f}")
    print(f"AUC Standard Deviation: {std_auc:.4f}")
    print(f"Minimum AUC recorded: {min(all_auc_scores):.4f}")
    print(f"Maximum AUC recorded: {max(all_auc_scores):.4f}")
    print("="*50)
    print(f"\n All 100 separate files were successfully generated inside:\n--> {output_subfolder}")

if __name__ == "__main__":
    main()