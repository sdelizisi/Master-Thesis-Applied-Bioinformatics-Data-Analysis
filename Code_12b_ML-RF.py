#!/usr/bin/env python3
# === ML Full Dataset Training with Random Forest and Resubstitution Evaluation ===
# This script loads the full transcriptomic matrix, trains a RandomForest classifier
# on 100% of the samples (N=40), evaluates resubstitution metrics without ROC/AUC,
# and exports the results to a text file.

# === Libraries ===
import random
from numpy import genfromtxt
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix

# === Seed Initialization ===
np.random.seed(1)
random.seed("Soultana")


# === Load Data Function ===
def load_data(dataset_path):
  dataset = genfromtxt(dataset_path, delimiter="\t", dtype="str")
  data = dataset[1:, 1:].astype("float")

  # Transpose matrix: scikit-learn expects samples as rows and genes as columns
  data = data.T

  print(f"Data Shape (Samples, Genes): {data.shape}")
  return data


# === Load Labels Function ===
def load_labels():
  # Strictly balanced cohort (N=40: 20 High, 20 Low Cortisol)
  labels = [1] * 20 + [0] * 20
  print(f"Labels loaded. Total samples: {len(labels)}")
  return np.array(labels)


# === Train Classifier Function ===
def train_classifier(X, y):
  classifier = RandomForestClassifier(criterion="gini", random_state=1)
  classifier.fit(X, y)
  return classifier


# === Evaluate Model Function (100% Resubstitution) ===
def evaluate_model(classifier, X, y):
  preds = classifier.predict(X)

  report = classification_report(y, preds, zero_division=0)
  print("\nClassification Report (100% Training Evaluation):")
  print(report)

  cm = confusion_matrix(y, preds)
  print(f"Confusion Matrix:\n{cm}")

  tn, fp, fn, tp = cm.ravel()
  specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
  sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
  accuracy = accuracy_score(y, preds)

  print(f"Accuracy:    {accuracy * 100:.2f}%")
  print(f"Sensitivity: {sensitivity * 100:.2f}%")
  print(f"Specificity: {specificity * 100:.2f}%")

  output_filename = "ML_Evaluation_Metrics_100_percent.txt"
  with open(output_filename, "w") as f:
    f.write("=== Random Forest Resubstitution Evaluation (100% Cohort) ===\n\n")
    f.write(f"Confusion Matrix:\n{cm}\n\n")
    f.write(f"Accuracy:    {accuracy * 100:.2f}%\n")
    f.write(f"Sensitivity: {sensitivity * 100:.2f}%\n")
    f.write(f"Specificity: {specificity * 100:.2f}%\n\n")
    f.write(f"Classification Report:\n{report}\n")

  print(f"Results successfully exported to '{output_filename}'.")


# === Main Function ===
def main():
  ranked_matrix_path = r"C:\Users\User\Desktop\Soultana\Msc Bioinformatics\Διπλωματική-Msc\07_Miscellaneous\Sfikakis_Yavropoulou\RNA-seq\files\Counted_Batch1_2_RawRead_PC_ranked.txt"

  data = load_data(ranked_matrix_path)
  labels = load_labels()

  # 100/0 strategy: train on the entire cohort
  classifier = train_classifier(data, labels)

  # Resubstitution performance evaluation on the training cohort
  evaluate_model(classifier, data, labels)


if __name__ == "__main__":
  main()