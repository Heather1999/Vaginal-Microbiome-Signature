import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import RocCurveDisplay, PrecisionRecallDisplay, roc_auc_score, average_precision_score
from sklearn.model_selection import GroupKFold, GridSearchCV
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression

# 1. Load and Prepare Data
meta_df = pd.read_csv("permanova/R_meta.tsv", sep="\t")
taxa_df = pd.read_csv("permanova/concated_clr.tsv", sep="\t", index_col=0)
# Add cohort metadata
meta_df["cohort"] = np.where(
    meta_df["Unnamed: 0"].astype(str).str.startswith("S"), "US", "DK"
)   
meta_df = meta_df.set_index("Unnamed: 0") 
# Join taxa and cohort info
X_df = taxa_df.join(meta_df[["cohort", "PatID"]])
y = meta_df["Archaea"].astype(int).values
print(X_df)
# --- STRIKE RIGOR STEP: SEPARATE US AND DK IMMEDIATELY ---
# We split the indices first
train_mask = X_df['cohort'] == 'US'
test_mask = X_df['cohort'] == 'DK'

X_train_full = X_df[train_mask].drop(columns=["PatID"])
y_train = y[train_mask]
groups_train = X_df[train_mask]["PatID"]

X_test_full = X_df[test_mask].drop(columns=["PatID"])
y_test = y[test_mask]

print(f"Training on US: {X_train_full.shape[0]} samples")
print(f"Testing on DK: {X_test_full.shape[0]} samples")

# 2. Setup Pipeline
taxa_cols = list(taxa_df.columns)
# Note: Since we are training ONLY on US, the 'cohort' column is constant (all "US")
# We remove 'cohort' from the model features to avoid errors in OneHotEncoding
preprocess = ColumnTransformer(
    transformers=[("taxa", StandardScaler(), taxa_cols)]
)

clf = LogisticRegression(
    penalty="elasticnet", solver="saga", max_iter=10000,
    class_weight="balanced", random_state=42
)

pipe = Pipeline([("prep", preprocess), ("clf", clf)])

# 3. GridSearchCV on US COHORT ONLY
cv = GroupKFold(n_splits=5)
param_grid = {
    "clf__C": [0.01, 0.1, 1, 10],
    "clf__l1_ratio": [0.1, 0.5, 0.9],
}

gs = GridSearchCV(
    pipe, param_grid=param_grid, scoring="roc_auc", 
    cv=cv, n_jobs=-1
)

# Fit only on US data
gs.fit(X_train_full, y_train, groups=groups_train)

print("\nBest params (tuned on US only):", gs.best_params_)
print("Internal US CV AUROC:", gs.best_score_)

# 4. Final External Validation on DK
best_model = gs.best_estimator_

# Predict on DK (The model has never seen this data or its distribution)
y_probs_dk = best_model.predict_proba(X_test_full)[:, 1]
external_auroc = roc_auc_score(y_test, y_probs_dk)
external_pr_auc = average_precision_score(y_test, y_probs_dk)

print(f"EXTERNAL VALIDATION (DK) AUROC: {external_auroc:.4f}")
print(f"EXTERNAL VALIDATION (DK) PR-AUC: {external_pr_auc:.4f}")

# 5. Visualization
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

RocCurveDisplay.from_predictions(y_test, y_probs_dk, name="US-trained Model", ax=axes[0])
axes[0].set_title("External Validation ROC (DK)")
axes[0].set_xlabel("False Positive Rate (Positive label: Archaea present)")
axes[0].set_ylabel("True Positive Rate (Positive label: Archaea present)")
axes[0].plot([0, 1], [0, 1], 'k--')

PrecisionRecallDisplay.from_predictions(y_test, y_probs_dk, name="US-trained Model", ax=axes[1])
axes[1].set_title("External Validation PR Curve (DK)")
axes[1].set_xlabel("Recall (Positive label: Archaea present)")
axes[1].set_ylabel("Precision (Positive label: Archaea present)")

plt.tight_layout()
plt.savefig("archaea_model_report.pdf", dpi=300, bbox_inches='tight', transparent=True)
plt.show()

# 6. Feature Selection (from the US-only model)
coefs = best_model.named_steps["clf"].coef_.ravel()
feature_names = best_model.named_steps["prep"].get_feature_names_out()

coef_df = pd.DataFrame({"feature": feature_names, "coefficient": coefs})
selected = coef_df[coef_df["coefficient"].abs() > 1e-6].sort_values("coefficient", ascending=False)
selected.to_csv("Strict_US_trained_Features_Archaea.tsv", sep='\t')