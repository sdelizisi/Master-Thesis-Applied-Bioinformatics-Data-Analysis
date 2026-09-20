#!/usr/bin/env python3
# === ML Predictions with Random Forest (Filtered Feature Importance > 0, 100/0 Strategy) ===
# This script filters the ranked dataset to keep ONLY genes with an Importance Score > 0,
# trains a RandomForest classifier on 100% of the data, and evaluates resubstitution performance.

# === Libraries ===
import os
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, roc_curve, auc
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import RocCurveDisplay

# === Seed Initialization ===
np.random.seed(1)
random.seed('Soultana')

# === Feature Selection: Get Active Genes Function ===
def get_active_genes(importance_file_path):
    print(f"Reading importance values from: {importance_file_path}")
    df_importance = pd.read_csv(importance_file_path, sep='\t')
    
    df_active = df_importance[df_importance['Importance_Score'] > 0]
    active_genes_list = df_active['Ensembl_ID'].tolist()
    
    print(f"--> Total genomic features with Gini Importance > 0: {len(active_genes_list)}")
    return active_genes_list

# === Load and Filter Data Function ===
def load_filtered_data(dataset_path, active_genes_list):
    print(f"Loading ranked expression matrix from: {dataset_path}")
    df = pd.read_csv(dataset_path, sep='\t', index_col=0)
    
    df_filtered = df.loc[df.index.isin(active_genes_list)]
    data = df_filtered.values.T
    
    print(f"--> Transposed Matrix Shape for ML (Samples, Filtered Genes): {data.shape}")
    return data

# === Load Labels Function ===
def load_labels():
    labels = [1] * 20 + [0] * 20
    print(f"Labels loaded. Total samples: {len(labels)}")
    return np.array(labels)

# === Train Classifier Function ===
def train_classifier(X_train, y_train):
    classifier = RandomForestClassifier(criterion='gini', random_state=1)
    classifier.fit(X_train, y_train)
    return classifier

# === Evaluate Model Function (Resubstitution / 100% Training) ===
def evaluate_model(classifier, X_train, y_train):
    
    probs = classifier.predict_proba(X_train)
    preds = probs[:, 1]
    
    fpr, tpr, _ = roc_curve(y_train, preds)
    roc_auc = auc(fpr, tpr)
    print(f"AUC Score (Resubstitution): {roc_auc:.4f}")

    # Ορισμός διαδρομής φακέλου και αυτόματη δημιουργία αν δεν υπάρχει
    target_dir = r"C:\Users\User\Desktop\Soultana\Msc Bioinformatics\Διπλωματική-Msc\07_Miscellaneous\Sfikakis_Yavropoulou\RNA-seq\files"
    os.makedirs(target_dir, exist_ok=True)

    # 1ο Διάγραμμα
    plt.figure()  
    RocCurveDisplay.from_predictions(
        y_train, 
        preds, 
        drop_intermediate=False, 
        name='Random Forest (100% Training)', 
        ax=plt.gca()
    )
    plt.title('ROC Curve (Resubstitution)')
    
    output_path_1 = os.path.join(target_dir, "AUC1_important.png")
    plt.savefig(output_path_1)
    plt.close()

    # Classification Report & Metrics
    mylist = (preds >= 0.5).astype(int)
    print("\nClassification Report (Resubstitution):")
    print(classification_report(y_train, mylist, zero_division=0))
    
    cm = confusion_matrix(y_train, mylist)
    if cm.size == 4:
        tn, fp, fn, tp = cm.ravel()
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
        print(f"Specificity: {specificity:.4f}")
        print(f"Sensitivity: {sensitivity:.4f}")
    
    accuracy = accuracy_score(y_train, mylist)
    print(f"Accuracy: {accuracy:.4f}")

    # 2ο Διάγραμμα
    plt.figure()  
    plt.title('Receiver Operating Characteristic (Resubstitution)')
    plt.plot(fpr, tpr, 'b', label=f'AUC = {roc_auc:.4f}')
    plt.legend(loc='lower right')
    plt.plot([0, 1], [0, 1], 'r--')
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.ylabel('True Positive Rate (Sensitivity)')
    plt.xlabel('False Positive Rate (1 - Specificity)')
    
    output_path_2 = os.path.join(target_dir, "AUC2_important.png")
    plt.savefig(output_path_2)
    plt.close() 
    
    print("ROC Curve plots successfully saved as 'AUC1_important.png' and 'AUC2_important.png' in your target directory.")

# === Main Function ===
def main():
    ranked_matrix_path = r"C:\Users\User\Desktop\Soultana\Msc Bioinformatics\Διπλωματική-Msc\07_Miscellaneous\Sfikakis_Yavropoulou\RNA-seq\files\Counted_Batch1_2_RawRead_PC_ranked.txt"
    importance_file_path = r"C:\Users\User\Desktop\Soultana\Msc Bioinformatics\Διπλωματική-Msc\07_Miscellaneous\Sfikakis_Yavropoulou\RNA-seq\files\Final_Table_Top_500_Genes.txt"
    
    active_genes = get_active_genes(importance_file_path)
    data = load_filtered_data(ranked_matrix_path, active_genes)
    labels = load_labels()
    
    X_train = data
    y_train = labels
    print(f"Training and evaluating model on 100% of the dataset (N = {len(y_train)} samples)...")
    
    classifier = train_classifier(X_train, y_train)
    evaluate_model(classifier, X_train, y_train)

if __name__ == "__main__":
    main()