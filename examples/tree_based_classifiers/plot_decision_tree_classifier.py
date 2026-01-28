"""
=========================================================
Decision Tree Classifier on Wine Dataset
=========================================================

This example demonstrates the usage of :class:`~sklearn.tree.DecisionTreeClassifier`
on the Wine dataset. We visualize the decision tree structure and evaluating
its performance.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import matplotlib.pyplot as plt

from sklearn.datasets import load_wine
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeClassifier, plot_tree

# %%
# Load Data
# ---------
# We generally load the wine dataset for this classification task.
data = load_wine()
X, y = data.data, data.target
feature_names = data.feature_names
class_names = data.target_names

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.3, random_state=42
)

# %%
# Train the Model
# ---------------
# We train a DecisionTreeClassifier with a maximum depth of 3 to keep the
# visualization manageable.
clf = DecisionTreeClassifier(max_depth=3, random_state=42)
clf.fit(X_train, y_train)

# %%
# Evaluate
# --------
# Check the accuracy on the test set.
y_pred = clf.predict(X_test)
accuracy = accuracy_score(y_test, y_pred)
print(f"Accuracy on test set: {accuracy:.2f}")

# %%
# Visualize the Tree
# ------------------
# We plot the tree structure to interpret how the model makes decisions.
plt.figure(figsize=(12, 8))
plot_tree(
    clf,
    feature_names=feature_names,
    class_names=class_names,
    filled=True,
    rounded=True,
    fontsize=10,
)
plt.title("Decision Tree trained on Wine Dataset")
plt.show()
