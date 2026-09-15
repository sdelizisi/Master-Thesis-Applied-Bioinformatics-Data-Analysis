#!/usr/bin/env python3
# === Nested Cross-Validation Script ===
# Applies 3x10 nested cross-validation for KNN, SVM, and Random Forest classifiers
# with hyperparameter tuning and evaluation.

# === Libraries ===
import random
import numpy as np
from numpy import genfromtxt
import pandas as pd
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# === Seed Initialization ===
np.random.seed(1)
random.seed('Soultana')

# === Load Data Function ===
def load_data(dataset_path):
    
    dataset = genfromtxt(dataset_path, delimiter='\t', dtype='str')
    data = dataset[1:, 1:].astype('float')
    
    # Transpose the data because scikit-learn requires samples as rows and genes as columns
    data = data.T 
    
    print(f"Data Shape (Samples, Genes): {data.shape}")
    return data

# === Load Labels Function ===
def load_labels():
    
    # 20 samples belonging to 'high' group (coded as 1)
    # 20 samples belonging to 'low' group (coded as 0)
    labels = [1] * 20 + [0] * 20
    
    print(f"Labels loaded. Total samples: {len(labels)}")
    return np.array(labels)

# === Nested Cross-Validation Function ===
def nested_cv(X, y, selected_classifiers):
   
    outer_cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=1)
    inner_cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=1)

    classifiers = {
        'KNN': (KNeighborsClassifier(), {
            'n_neighbors': list(range(1, 15)), # Adjusted max range for n=40 samples
            'weights': ['uniform', 'distance'],
            'metric': ['manhattan', 'euclidean', 'chebyshev', 'minkowski']
        }),
        'SVM': (SVC(probability=True), {
            'C': [0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000],
            'kernel': ['linear', 'rbf'],
            'gamma': ['scale', 0.001, 0.0001]
        }),
        'RF': (RandomForestClassifier(random_state=1), {
            'criterion': ['gini', 'entropy']
        })
    }

    
    active_classifiers = {name: classifiers[name] for name in selected_classifiers if name in classifiers}

    for clf_name, (clf, param_grid) in active_classifiers.items():
        print(f"=== Evaluating Classifier: {clf_name} ===")
        outer_accuracies = []

        for outer_fold, (train_idx, test_idx) in enumerate(outer_cv.split(X, y), start=1):
            print(f"\n--- Outer Fold {outer_fold} ---")
            X_train, X_test = X[train_idx], X[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]

            # Inner CV with grid search
            grid_search = GridSearchCV(clf, param_grid, cv=inner_cv, scoring='accuracy', n_jobs=-1)
            grid_search.fit(X_train, y_train)

            best_model = grid_search.best_estimator_
            predictions = best_model.predict(X_test)

            acc = accuracy_score(y_test, predictions)
            outer_accuracies.append(acc)

            print(f"Best Params: {grid_search.best_params_}")
            print(f"Accuracy (Outer Fold {outer_fold}): {acc:.4f}")
            print("Classification Report:")
            print(classification_report(y_test, predictions, zero_division=0))

        print(f"\n>>> FINAL OVERALL RESULTS FOR {clf_name} <<<")
        print(f"Mean Accuracy for {clf_name}: {np.mean(outer_accuracies):.4f}")
        print(f"Standard Deviation: {np.std(outer_accuracies):.4f}")

# === Main Function ===
def main():
    # Load the data and labels using your precise file path
    ranked_matrix_path = r"C:\Users\User\Desktop\Soultana\Msc Bioinformatics\Διπλωματική-Msc\07_Miscellaneous\Sfikakis_Yavropoulou\RNA-seq\files\Counted_Batch1_2_RawRead_PC_ranked.txt"
    
    data = load_data(ranked_matrix_path)
    labels = load_labels()

    
    run_knn = False    
    run_svm = False  
    run_rf  = True  

    selected_classifiers = []
    if run_knn: selected_classifiers.append('KNN')
    if run_svm: selected_classifiers.append('SVM')
    if run_rf:  selected_classifiers.append('RF')

    if len(selected_classifiers) == 0:
        print("\n[!] Error! Choose an algorithm!")
        return

    # Perform Nested Cross-Validation
    nested_cv(data, labels, selected_classifiers)

# === Run the Script ===
if __name__ == "__main__":
    main()










