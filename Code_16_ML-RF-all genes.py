#!/usr/bin/env python3
# === ML Predictions with Random Forest (ALL Genes, 100/0 Strategy) ===
# This script trains a RandomForest classifier on 100% of the data using ALL genomic features, 
# and evaluates resubstitution performance.

# === Libraries ===
import random
from numpy import genfromtxt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, roc_curve, auc
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import RocCurveDisplay

# === Seed Initialization ===
# Set seeds for reproducibility
np.random.seed(1)
random.seed('Soultana')

# === Load All Data Function ===
def load_all_data(dataset_path):
    print(f"Loading full ranked expression matrix from: {dataset_path}")
    # Load with index_col=0 to cleanly handle Ensembl IDs as row names
    df = pd.read_csv(dataset_path, sep='\t', index_col=0)
    
    # Convert values to a float matrix and transpose (Samples as rows, Genes as columns)
    data = df.values.T
    
    print(f"--> Transposed Matrix Shape for ML (Samples, ALL Genes): {data.shape}")
    return data

# === Load Labels Function ===
def load_labels():
    # Adjusted strictly for n=40 samples to match the DESeq2 metadata design
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
    
    # Get probabilities for ROC curve on training data
    probs = classifier.predict_proba(X_train)
    preds = probs[:, 1]
    
    # Calculate ROC AUC
    fpr, tpr, _ = roc_curve(y_train, preds)
    roc_auc = auc(fpr, tpr)
    print(f"AUC Score (Resubstitution - All Genes): {roc_auc:.4f}")

    # ROC curve plot using sklearn tool 
    plt.figure()  
    RocCurveDisplay.from_predictions(
        y_train, 
        preds, 
        drop_intermediate=False, 
        name='Random Forest (All Genes)', 
        ax=plt.gca()
    )
    plt.title('ROC Curve (All Genes - Resubstitution)')
    
    # Save 1st plot with new name and directory
    output_path_1 = r"C:\Users\User\Desktop\Soultana\Msc Bioinformatics\Διπλωματική-Msc\07_Miscellaneous\Sfikakis_Yavropoulou\RNA-seq\files\AUC1_all_genes.png"
    plt.savefig(output_path_1)
    plt.close()

    # Calculate confusion matrix
    mylist = (preds >= 0.5).astype(int)  # Thresholding prediction at 0.5
    print("\nClassification Report (Resubstitution - All Genes):")
    print(classification_report(y_train, mylist, zero_division=0))
    
    # Confusion Matrix Parsing
    cm = confusion_matrix(y_train, mylist)
    if cm.size == 4:
        tn, fp, fn, tp = cm.ravel()
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
        print(f"Specificity: {specificity:.4f}")
        print(f"Sensitivity: {sensitivity:.4f}")
    
    accuracy = accuracy_score(y_train, mylist)
    print(f"Accuracy: {accuracy:.4f}")

    # Generate and save the custom ROC curve plot
    plt.figure()  
    plt.title('Receiver Operating Characteristic (All Genes)')
    plt.plot(fpr, tpr, 'b', label=f'AUC = {roc_auc:.4f}')
    plt.legend(loc='lower right')
    plt.plot([0, 1], [0, 1], 'r--')
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.ylabel('True Positive Rate (Sensitivity)')
    plt.xlabel('False Positive Rate (1 - Specificity)')
    
    # Save 2nd plot with new name and directory
    output_path_2 = r"C:\Users\User\Desktop\Soultana\Msc Bioinformatics\Διπλωματική-Msc\07_Miscellaneous\Sfikakis_Yavropoulou\RNA-seq\files\AUC2_all_genes.png"
    plt.savefig(output_path_2)
    plt.close() 
    
    print("ROC Curve plots successfully saved as 'AUC1_all_genes.png' and 'AUC2_all_genes.png' in your target directory.")

# === Main Function ===
def main():
    
    # Absolute path to full ranked dataset
    ranked_matrix_path = r"C:\Users\User\Desktop\Soultana\Msc Bioinformatics\Διπλωματική-Msc\07_Miscellaneous\Sfikakis_Yavropoulou\RNA-seq\files\Counted_Batch1_2_RawRead_PC_ranked.txt"
    
    # 1. Load full expression matrix (All genes)
    data = load_all_data(ranked_matrix_path)
    labels = load_labels()
    
    # 2. Full training strategy (100% training, 0% test split / Resubstitution)
    X_train = data
    y_train = labels
    print(f"Training and evaluating model on 100% of the dataset with ALL genes (N = {len(y_train)} samples)...")
    
    # 3. Train the classifier on the full feature space
    classifier = train_classifier(X_train, y_train)

    # 4. Evaluate the model on the training set
    evaluate_model(classifier, X_train, y_train)

# Run the main function
if __name__ == "__main__":
    main()