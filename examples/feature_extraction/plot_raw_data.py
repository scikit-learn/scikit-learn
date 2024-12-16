# %%
"""
============================
Feature Extraction Examples
============================

This example demonstrates how to use the `DictVectorizer` and `FeatureHasher` from the
`sklearn.feature_extraction` module to transform a list of dictionaries
into a feature matrix.

"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause
# %%
from sklearn.feature_extraction import DictVectorizer, FeatureHasher

# Sample data: list of dictionaries
data = [
    {"feature1": 1, "feature2": 2},
    {"feature1": 3, "feature2": 4},
    {"feature1": 5, "feature2": 6},
]

# Initialize the DictVectorizer
vec = DictVectorizer()

# Transform the data into a feature matrix using DictVectorizer
feature_matrix = vec.fit_transform(data).toarray()

# Print the feature names
print("Feature names:", vec.get_feature_names_out())

# Print the feature matrix
print("Feature matrix:\n", feature_matrix)

# Initialize the FeatureHasher
hasher = FeatureHasher(n_features=10, input_type="dict")

# Transform the data into a hashed feature matrix using FeatureHasher
hashed_matrix = hasher.transform(data).toarray()

# Print the hashed feature matrix
print("Hashed feature matrix:\n", hashed_matrix)
