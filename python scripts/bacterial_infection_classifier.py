"""
Build classifiers to predict bacterial vs non-bacterial infection from nasal swab gene expression data.

Uses a Bagged Support Vector Machine (BagSVM) with nested cross-validation:
- Outer CV: Model evaluation
- Inner CV: Hyperparameter tuning (GridSearchCV) and feature selection (RFECV)

Input files:
    - Gene counts CSV (samples x genes)
    - Metadata CSV with target labels (bacterial vs non-bacterial)
    - Differentially expressed genes CSV

Output files:
    - Trained classifier dump files (.joblib)
    - Selected feature/predictor lists
    - Performance summary table with AUC-ROC scores
"""

import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC
from sklearn.ensemble import BaggingClassifier
from sklearn.feature_selection import RFECV
from joblib import dump, load

# Configuration
RESULTS_PATH = "path/to/folder"
OUTPUT_PREFIX = "bacterial_infection_classifier"
MODE = "create"  # "create" to train new models or "load" to load existing
NUM_CV = 5  # Number of cross-validation folds
TEST_PROPORTION = 0.25  # Proportion of data held out for final testing
NUM_OUTER_CV = 10  # Number of outer CV iterations

# Load data
gene_counts = pd.read_csv("genecounts.csv",index_col=0)

metadata = pd.read_csv("metadata.csv",dtype={'ID': np.object_}).set_index("ID")

deg_results = pd.read_csv("deg.csv", index_col=0)

# Align metadata with gene count samples
metadata = metadata.loc[gene_counts.columns, :]

# Get gene names and indices of differentially expressed genes
all_gene_names = gene_counts.index.values
deg_indices = [all_gene_names.tolist().index(gene) for gene in deg_results.index.values]

# Set up the machine learning pipeline
# RFECV: Recursive feature elimination with cross-validation
# Drops 0.1 -> 10% of features per iteration
# Always include at least 2 predictors
rfecv = RFECV(
    estimator=BaggingClassifier(estimator=LinearSVC(max_iter=10000), random_state=123),
    cv=NUM_CV,
    n_jobs=1,  # Set to 1 to avoid nested parallelism issues
    scoring='roc_auc',
    step=0.1,  # Drop 10% of features per iteration
    min_features_to_select=2,
    verbose=True
)

# Complete pipeline: scale → feature selection → classification
pipeline = Pipeline([
    ('scale', StandardScaler()),
    ('rfe', rfecv),
    ('classifier', BaggingClassifier(estimator=LinearSVC(max_iter=10000)))
])

# Hyperparameter grid for tuning
param_grid = {
    'classifier__random_state': [123],
    'classifier__base_estimator__C': [0.01, 0.1, 1, 10],  # Regularization strength
    'classifier__max_features': [0.1, 0.5, 0.8],  # Proportion of features per estimator
    'classifier__n_estimators': [100, 1000]  # Number of base estimators
}

# Grid search with cross-validation
grid_search = GridSearchCV(
    estimator=pipeline,
    cv=NUM_CV,
    n_jobs=1,
    scoring='roc_auc',
    param_grid=param_grid,
    verbose=True
)

# Prepare data
# Transpose so rows are samples, columns are genes
gene_counts = gene_counts.T
target_labels = metadata.bv_both_none  # Bacterial vs non-bacterial labels

# Create held-out test set (stratified split)
X_train_full, X_test_full, y_train_full, y_test_full = train_test_split(
    gene_counts,
    target_labels,
    test_size=TEST_PROPORTION,
    random_state=123,
    stratify=target_labels
)

# Storage for cross-validation results
cv_results_train = []
cv_results_test = []
cv_results_test_full = []

# Outer cross-validation loop
for cv_iter in range(NUM_OUTER_CV):
    cv_id = cv_iter + 1
    print(f"\n{'='*50}\nCV Iteration {cv_id}/{NUM_OUTER_CV}\n{'='*50}")
    
    # Split training data for this CV iteration
    X_train, X_test, y_train, y_test = train_test_split(
        X_train_full,
        y_train_full,
        test_size=TEST_PROPORTION,
        stratify=y_train_full
    )
    
    # Train or load model
    model_path = f"{RESULTS_PATH}/{OUTPUT_PREFIX}_cv{cv_id}_model.joblib"
    
    if MODE == "create":
        grid_search.fit(X_train, y_train)
        dump(grid_search, model_path)
    else:
        grid_search = load(model_path)
    
    # Extract selected features
    selected_features = all_gene_names[deg_indices][
        grid_search.best_estimator_.named_steps["rfe"].support_
    ]
    
    # Save selected features
    pd.DataFrame(selected_features).to_csv(
        f"{RESULTS_PATH}/{OUTPUT_PREFIX}_cv{cv_id}_features.csv",
        header=False,
        index=False
    )
    
    print(f"Best parameters: {grid_search.best_params_}")
    print(f"Selected features ({len(selected_features)}): {selected_features[:5]}...")
    
    # Evaluate on training set
    train_probs = grid_search.predict_proba(X_train)[:, 1]
    train_auc = roc_auc_score(y_train, train_probs)
    cv_results_train.append({
        "labels": y_train,
        "probs": train_probs,
        "n_features": len(selected_features),
        "roc_auc": train_auc
    })
    print(f"Train AUC: {train_auc:.3f}")
    
    # Evaluate on test set (from this CV split)
    test_probs = grid_search.predict_proba(X_test)[:, 1]
    test_auc = roc_auc_score(y_test, test_probs)
    cv_results_test.append({
        "labels": y_test,
        "probs": test_probs,
        "roc_auc": test_auc
    })
    print(f"Test AUC: {test_auc:.3f}")
    
    # Evaluate on held-out test set
    test_full_probs = grid_search.predict_proba(X_test_full)[:, 1]
    test_full_auc = roc_auc_score(y_test_full, test_full_probs)
    cv_results_test_full.append({
        "labels": y_test_full,
        "probs": test_full_probs,
        "roc_auc": test_full_auc
    })
    print(f"Held-out test AUC: {test_full_auc:.3f}")

# Train final model on full training set
print(f"\n{'='*50}\nTraining final model on full training set\n{'='*50}")
final_model_path = f"{RESULTS_PATH}/{OUTPUT_PREFIX}_final_model.joblib"

if MODE == "create":
    grid_search.fit(X_train_full, y_train_full)
    dump(grid_search, final_model_path)
else:
    grid_search = load(final_model_path)

# Extract and save final features
final_features = all_gene_names[deg_indices][
    grid_search.best_estimator_.named_steps["rfe"].support_
]
pd.DataFrame(final_features).to_csv(
    f"{RESULTS_PATH}/{OUTPUT_PREFIX}_final_features.csv",
    header=False,
    index=False
)

# Evaluate final model
final_train_probs = grid_search.predict_proba(X_train_full)[:, 1]
final_train_auc = roc_auc_score(y_train_full, final_train_probs)

final_test_probs = grid_search.predict_proba(X_test_full)[:, 1]
final_test_auc = roc_auc_score(y_test_full, final_test_probs)

print(f"Final model features: {len(final_features)}")
print(f"Final train AUC: {final_train_auc:.3f}")
print(f"Final held-out test AUC: {final_test_auc:.3f}")

# Create performance summary table
summary_data = {
    'cv_id': list(range(1, NUM_OUTER_CV + 1)),
    'n_features': [x["n_features"] for x in cv_results_train],
    'train_auc': [round(x["roc_auc"], 3) for x in cv_results_train],
    'test_auc': [round(x["roc_auc"], 3) for x in cv_results_test],
    'held_out_auc': [round(x["roc_auc"], 3) for x in cv_results_test_full]
}

# Add mean and standard deviation row
train_aucs = summary_data['train_auc']
test_aucs = summary_data['test_auc']
held_out_aucs = summary_data['held_out_auc']

summary_data['cv_id'].append("Mean (SD)")
summary_data['n_features'].append("")
summary_data['train_auc'].append(f"{np.mean(train_aucs):.3f} ({np.std(train_aucs):.3f})")
summary_data['test_auc'].append(f"{np.mean(test_aucs):.3f} ({np.std(test_aucs):.3f})")
summary_data['held_out_auc'].append(f"{np.mean(held_out_aucs):.3f} ({np.std(held_out_aucs):.3f})")

# Add final model row
summary_data['cv_id'].append("Final")
summary_data['n_features'].append(len(final_features))
summary_data['train_auc'].append(round(final_train_auc, 3))
summary_data['test_auc'].append("N/A")
summary_data['held_out_auc'].append(round(final_test_auc, 3))

# Save summary
summary_df = pd.DataFrame(summary_data)
summary_df.to_csv(
    f"{RESULTS_PATH}/{OUTPUT_PREFIX}_performance_summary.csv",
    index=False
)

print("\n" + "="*50)
print("Performance Summary:")
print("="*50)
print(summary_df.to_string(index=False))
