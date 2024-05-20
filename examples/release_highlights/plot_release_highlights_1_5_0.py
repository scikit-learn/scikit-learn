# ruff: noqa
"""
=======================================
Release Highlights for scikit-learn 1.5
=======================================

.. currentmodule:: sklearn

We are pleased to announce the release of scikit-learn 1.5! Many bug fixes
and improvements were added, as well as some key new features. Below we
detail the highlights of this release. **For an exhaustive list of
all the changes**, please refer to the :ref:`release notes <release_notes_1_5>`.

To install the latest version (with pip)::

    pip install --upgrade scikit-learn

or with conda::

    conda install -c conda-forge scikit-learn

"""

# %%
# FixedThresholdClassifier: Setting the decision threshold of a binary classifier
# -------------------------------------------------------------------------------
# All binary classifiers of scikit-learn use a fixed decision threshold of 0.5 to
# convert probability estimates (i.e. output of `predict_proba`) into class
# predictions. However, 0.5 is almost never the desired threshold for a given problem.
# :class:`~model_selection.FixedThresholdClassifier` allows to wrap any binary
# classifier and set a custom decision threshold.
from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix

X, y = make_classification(n_samples=1_000, weights=[0.9, 0.1], random_state=0)
classifier = LogisticRegression(random_state=0).fit(X, y)

print("confusion matrix:\n", confusion_matrix(y, classifier.predict(X)))

# %%
# Lowering the threshold, i.e. allowing more samples to be classified as the positive
# class, increases the number of true positives at the cost of more false positives
# (as is well known from the concavity of the ROC curve).
from sklearn.model_selection import FixedThresholdClassifier

wrapped_classifier = FixedThresholdClassifier(classifier, threshold=0.1).fit(X, y)

print("confusion matrix:\n", confusion_matrix(y, wrapped_classifier.predict(X)))

# %%
# TunedThresholdClassifierCV: Tuning the decision threshold of a binary classifier
# --------------------------------------------------------------------------------
# The decision threshold of a binary classifier can be tuned to optimize a given
# metric, using :class:`~model_selection.TunedThresholdClassifierCV`.
from sklearn.metrics import balanced_accuracy_score

# Due to the class imbalance, the balanced accuracy is not optimal for the default
# threshold. The classifier tends to over predict the majority class.
print(f"balanced accuracy: {balanced_accuracy_score(y, classifier.predict(X)):.2f}")

# %%
# Tuning the threshold to optimize the balanced accuracy gives a smaller threshold
# that allows more samples to be classified as the positive class.
from sklearn.model_selection import TunedThresholdClassifierCV

tuned_classifier = TunedThresholdClassifierCV(
    classifier, cv=5, scoring="balanced_accuracy"
).fit(X, y)

print(f"new threshold: {tuned_classifier.best_threshold_:.4f}")
print(
    f"balanced accuracy: {balanced_accuracy_score(y, tuned_classifier.predict(X)):.2f}"
)

# %%
# :class:`~model_selection.TunedThresholdClassifierCV` also benefits from the
# metadata routing support (:ref:`Metadata Routing User Guide<metadata_routing>`)
# allowing to optimze complex business metrics, detailed
# in :ref:`Post-tuning the decision threshold for cost-sensitive learning
# <sphx_glr_auto_examples_model_selection_plot_cost_sensitive_learning.py>`.

# %%
# Performance improvements in PCA
# -------------------------------
# :class:`~decomposition.PCA` has a new solver, "covariance_eigh", which is faster
# and more memory efficient than the other solvers for datasets with a large number
# of samples and a small number of features.
from sklearn.datasets import make_low_rank_matrix
from sklearn.decomposition import PCA

X = make_low_rank_matrix(
    n_samples=10_000, n_features=100, tail_strength=0.1, random_state=0
)

pca = PCA(n_components=10).fit(X)

print(f"explained variance: {pca.explained_variance_ratio_.sum():.2f}")

# %%
# The "full" solver has also been improved to use less memory and allows to
# transform faster. The "auto" option for the solver takes advantage of the
# new solver and is now able to select an appropriate solver for sparse
# datasets.
from scipy.sparse import random

X = random(10000, 100, format="csr", random_state=0)

pca = PCA(n_components=10, svd_solver="auto").fit(X)

# %%
# ColumnTransformer is subscriptable
# ----------------------------------
# The transformers of a :class:`~compose.ColumnTransformer` can now be directly
# accessed using indexing by name.
import numpy as np
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler, OneHotEncoder

X = np.array([[0, 1, 2], [3, 4, 5]])
column_transformer = ColumnTransformer(
    [("std_scaler", StandardScaler(), [0]), ("one_hot", OneHotEncoder(), [1, 2])]
)

column_transformer.fit(X)

print(column_transformer["std_scaler"])
print(column_transformer["one_hot"])

# %%
# Custom imputation strategies for the SimpleImputer
# --------------------------------------------------
# :class:`~impute.SimpleImputer` now supports custom strategies for imputation,
# using a callable that computes a scalar value from the non missing values of
# a column vector.
from sklearn.impute import SimpleImputer

X = np.array(
    [
        [-1.1, 1.1, 1.1],
        [3.9, -1.2, np.nan],
        [np.nan, 1.3, np.nan],
        [-0.1, -1.4, -1.4],
        [-4.9, 1.5, -1.5],
        [np.nan, 1.6, 1.6],
    ]
)


def smallest_abs(arr):
    """Return the smallest absolute value of a 1D array."""
    return np.min(np.abs(arr))


imputer = SimpleImputer(strategy=smallest_abs)

imputer.fit_transform(X)

# %%
# Pairwise distances with non-numeric arrays
# ------------------------------------------
# :func:`~metrics.pairwise_distances` can now compute distances between
# non-numeric arrays using a callable metric.
from sklearn.metrics import pairwise_distances

X = ["cat", "dog"]
Y = ["cat", "fox"]


def levenshtein_distance(x, y):
    """Return the Levenshtein distance between two strings."""
    if x == "" or y == "":
        return max(len(x), len(y))
    if x[0] == y[0]:
        return levenshtein_distance(x[1:], y[1:])
    return 1 + min(
        levenshtein_distance(x[1:], y),
        levenshtein_distance(x, y[1:]),
        levenshtein_distance(x[1:], y[1:]),
    )


pairwise_distances(X, Y, metric=levenshtein_distance)
