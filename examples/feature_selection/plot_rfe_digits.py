"""
=============================
Recursive feature elimination
=============================

This example demonstrates how :class:`~sklearn.feature_selection.RFE` can be used
to determine the importance of individual pixels when classifying handwritten digits.
RFE is a method that recursively removes the least significant features and retrains
the model, allowing us to rank features by their importance.

.. note::

    See also :ref:`sphx_glr_auto_examples_feature_selection_plot_rfe_with_cross_validation.py`

"""  # noqa: E501

# %%
# Dataset
# -------
#
# We will start by loading the handwritten digits dataset. This dataset consists of 8x8
# pixel images of handwritten digits. Each pixel will be treated as a feature and we
# aim to determine which pixels are most relevant for the digit classification task.

# %%
import matplotlib.pyplot as plt

from sklearn.datasets import load_digits

# Load the digits dataset
digits = load_digits()
X = digits.images.reshape((len(digits.images), -1))
y = digits.target

# Display the first digit
plt.imshow(digits.images[0], cmap="gray")
plt.title(f"Label: {digits.target[0]}")
plt.axis("off")
plt.show()

# %%
# Splitting the dataset for evaluation
# ------------------------------------
#
# To assess the benefits of feature selection with
# :class:`~sklearn.feature_selection.RFE`, we need a training set for selecting
# features and training our model, and a test set for evaluation.
# We'll allocate 70% of the data for training and 30% for testing.

# %%
from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.3, random_state=42
)

# %%
# Benchmarking SVM without Feature Selection
# ------------------------------------------
#
# Before applying :class:`~sklearn.feature_selection.RFE`, let's benchmark the
# performance of a :class:`~sklearn.svm.SVC` using all features. This will give us
# a baseline accuracy to compare against.

# %%
from sklearn.metrics import accuracy_score
from sklearn.svm import SVC

svc = SVC(kernel="linear", C=1)
svc.fit(X_train, y_train)
y_pred = svc.predict(X_test)
accuracy_all_features = accuracy_score(y_test, y_pred)

print(f"Accuracy using all {X_train.shape[1]} features: {accuracy_all_features:.4f}")

# %%
# Feature Selection with RFE
# --------------------------
#
# Now, we'll employ :class:`~sklearn.feature_selection.RFE` to select a subset of
# the most discriminative features. The goal is to determine if a reduced set of
# important features can either maintain or even improve the classifier's performance.

# %%
from sklearn.feature_selection import RFE

# Arbitrarily chosen; can be adjusted based on domain knowledge or iterative testing
num_features_to_select = 40
rfe = RFE(estimator=svc, n_features_to_select=num_features_to_select, step=10)
rfe.fit(X_train, y_train)

# %%
# Evaluating SVM on Selected Features
# -----------------------------------
#
# With the top features selected by :class:`~sklearn.feature_selection.RFE`, let's
# train a new :class:`~sklearn.svm.SVC` and assess its performance. The idea is to
# observe if there's any significant change in accuracy, ideally aiming for improvement.

# %%
X_train_rfe = rfe.transform(X_train)
X_test_rfe = rfe.transform(X_test)

svc_rfe = SVC(kernel="linear", C=1)
svc_rfe.fit(X_train_rfe, y_train)
y_pred_rfe = svc_rfe.predict(X_test_rfe)
accuracy_selected_features = accuracy_score(y_test, y_pred_rfe)

print(
    f"Accuracy using {num_features_to_select} selected features:"
    f" {accuracy_selected_features:.4f}"
)

# %%
# Visualizing Feature Importance after RFE
# ----------------------------------------
#
# :class:`~sklearn.feature_selection.RFE` provides a ranking of the features based on
# their importance. We can visualize this ranking to gain insights into which pixels
# (or features) are deemed most significant by :class:`~sklearn.feature_selection.RFE`
# in the digit classification task.

# %%
ranking = rfe.ranking_.reshape(digits.images[0].shape)
plt.matshow(ranking, cmap=plt.cm.Blues)
plt.colorbar()
plt.title("Ranking of pixels with RFE")
plt.show()

# %%
# Feature Selection Impact on Model Accuracy
# ---------------------------------------------------
#
# To understand the relationship between the number of features selected and model
# performance, let's train the :class:`~sklearn.svm.SVC` on various subsets of
# features ranked by :class:`~sklearn.feature_selection.RFE`. We'll then plot the
# accuracy of the model as a function of the number of features used. This will help
# us visualize any trade-offs between feature selection and model accuracy.

# %%
import numpy as np

# Split the dataset
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.3, random_state=42
)

# Train with RFE to get the rankings (as done earlier in the code)
svc = SVC(kernel="linear", C=1)
rfe = RFE(estimator=svc, n_features_to_select=1, step=1)
rfe.fit(X_train, y_train)
ranking = rfe.ranking_

# Store accuracies
num_features_list = [
    1,
    5,
    10,
    20,
    30,
    40,
    50,
    64,
]  # Adjust the step for finer granularity
accuracies = []

for num_features in num_features_list:
    # Select top 'num_features' important features
    top_features_idx = np.where(ranking <= num_features)[0]
    X_train_selected = X_train[:, top_features_idx]
    X_test_selected = X_test[:, top_features_idx]

    # Train SVM and get accuracy
    svc_selected = SVC(kernel="linear", C=1)
    svc_selected.fit(X_train_selected, y_train)
    y_pred = svc_selected.predict(X_test_selected)
    accuracy = accuracy_score(y_test, y_pred)
    accuracies.append(accuracy)

# Plot the accuracies
plt.plot(num_features_list, accuracies, marker="o", linestyle="-")
plt.xlabel("Number of Selected Features")
plt.ylabel("Accuracy")
plt.title("Feature Selection Impact on Model Accuracy")
plt.grid(True)
plt.show()
