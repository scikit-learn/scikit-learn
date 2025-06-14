"""
==========================
Metadata Routing Tree View
==========================

.. currentmodule:: sklearn

This example demonstrates the improved tree-based visualization for metadata routing
that makes it easier to understand how parameters are routed through a complex
pipeline of estimators.

The tree view reduces verbosity by showing each component only once and displaying
parameters and method mappings in a more structured way.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%

from sklearn import set_config
from sklearn.compose import ColumnTransformer
from sklearn.feature_selection import SelectPercentile, chi2
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold, RandomizedSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.utils._metadata_routing_tree_visualise import visualise_routing_tree

# Import both visualization methods for comparison
from sklearn.utils._metadata_routing_visualise import (
    visualise_routing as original_visualise,
)
from sklearn.utils.metadata_routing import get_routing_for_object

# Enable metadata routing
set_config(enable_metadata_routing=True)

# Create a complex pipeline similar to the one in the original example
numeric_features = ["age", "fare"]
numeric_transformer = Pipeline(
    steps=[
        ("imputer", SimpleImputer(strategy="median")),
        (
            "scaler",
            StandardScaler()
            .set_fit_request(sample_weight=True)
            .set_transform_request(copy=True),
        ),
    ]
)

categorical_features = ["embarked", "sex", "pclass"]
categorical_transformer = Pipeline(
    steps=[
        ("encoder", OneHotEncoder(handle_unknown="ignore")),
        ("selector", SelectPercentile(chi2, percentile=50)),
    ]
)
preprocessor = ColumnTransformer(
    transformers=[
        ("num", numeric_transformer, numeric_features),
        ("cat", categorical_transformer, categorical_features),
    ]
)

# Append classifier to preprocessing pipeline
clf = Pipeline(
    steps=[
        ("preprocessor", preprocessor),
        ("classifier", LogisticRegression().set_fit_request(sample_weight=True)),
    ]
)

param_grid = {
    "preprocessor__num__imputer__strategy": ["mean", "median"],
    "preprocessor__cat__selector__percentile": [10, 30, 50, 70],
    "classifier__C": [0.1, 1.0, 10, 100],
}

search_cv = RandomizedSearchCV(clf, param_grid, cv=GroupKFold(), random_state=0)

# Get the routing information
routing_info = get_routing_for_object(search_cv)

# Compare the two visualization methods
print("=" * 80)
print("ORIGINAL VERBOSE VISUALIZATION")
print("=" * 80)
# Use the original visualization function
original_visualise(routing_info)

print("\n" + "=" * 80)
print("NEW TREE-BASED VISUALIZATION")
print("=" * 80)
print(visualise_routing_tree(routing_info))

print("\n" + "=" * 80)
print("PARAMETER-SPECIFIC VISUALIZATION: sample_weight")
print("=" * 80)
print(visualise_routing_tree(routing_info, param="sample_weight"))

print("\n" + "=" * 80)
print("PARAMETER-SPECIFIC VISUALIZATION: copy")
print("=" * 80)
print(visualise_routing_tree(routing_info, param="copy"))

print("\n" + "=" * 80)
print("PARAMETER-SPECIFIC VISUALIZATION: groups")
print("=" * 80)
print(visualise_routing_tree(routing_info, param="groups"))
