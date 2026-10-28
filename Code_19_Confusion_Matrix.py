import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, roc_auc_score

# File paths
FILE_100 = r'C:\Users\User\Desktop\Soultana\Msc Bioinformatics\Διπλωματική-Msc\07_Miscellaneous\Sfikakis_Yavropoulou\RNA-seq\files\Counted_Batch1_2_RawRead_PC_ranked.txt'
FILE_90 = os.path.join(os.path.dirname(FILE_100), 'Counted_Batch1_2_RawRead_PC_ranked_Symbols.txt')

y = np.array([0]*20 + [1]*20)  # 0 = Low Cortisol, 1 = High Cortisol
names = ["Low Cortisol", "High Cortisol"]

# --- 1. 100/0 Resubstitution (full set of 40 samples) ---
X = pd.read_csv(FILE_100, sep="\t", index_col=0).T
rf_100 = RandomForestClassifier(criterion="gini", random_state=1, n_jobs=-1)
rf_100.fit(X, y)
cm_100 = confusion_matrix(y, rf_100.predict(X), labels=[0, 1])

# --- 2. 90/10: sum of 100 runs (objective estimate of generalization) ---
X = pd.read_csv(FILE_90, sep="\t", index_col=0).T
cm_90 = np.zeros((2, 2), dtype=int)
aucs = []

for i in range(1, 101):
    X_tr, X_te, y_tr, y_te = train_test_split(X, y, test_size=0.1, stratify=y, random_state=i)
    rf_90 = RandomForestClassifier(random_state=i, n_jobs=-1)
    rf_90.fit(X_tr, y_tr)
    cm_90 += confusion_matrix(y_te, rf_90.predict(X_te), labels=[0, 1])
    aucs.append(roc_auc_score(y_te, rf_90.predict_proba(X_te)[:, 1]))

print("Mean AUC:", round(np.mean(aucs), 4), "| SD:", round(np.std(aucs, ddof=1), 4))
print("90/10 aggregated confusion matrix (100 runs):\n", cm_90)

# --- 3. Visualization: TN/FP/FN/TP label + % per true class + absolute count ---
cell_labels = np.array([["TN", "FP"],
                        ["FN", "TP"]])

fig, ax = plt.subplots(1, 2, figsize=(12, 5))
titles = ["100/0 Resubstitution (N = 40)", "90/10, sum of 100 runs (N = 400 predictions)"]
for a, cm, title_text in zip(ax, [cm_100, cm_90], titles):
    cm_pct = cm / cm.sum(axis=1, keepdims=True) * 100
    disp = ConfusionMatrixDisplay(confusion_matrix=cm_pct, display_labels=names)
    disp.plot(ax=a, cmap="Blues", colorbar=False)
    for r in range(2):
        for c in range(2):
            t = disp.text_[r, c]
            t.set_text(f"{cell_labels[r, c]}\n{cm_pct[r, c]:.1f}%\n(n = {cm[r, c]})")
            t.set_fontsize(12)
            t.set_fontweight('bold')
    a.set_title(title_text, fontsize=13, fontweight='bold', pad=12)

plt.tight_layout()
plt.savefig("confusion_matrices_final.png", dpi=300)
plt.show()