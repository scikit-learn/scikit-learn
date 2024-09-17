"""
===================================================
Recursive feature elimination with cross-validation
===================================================

A Recursive Feature Elimination (RFE) example with automatic tuning of the
number of features selected with cross-validation.

"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Data generation
# ---------------
#
# We build a classification task using 3 informative features. The introduction
# of 2 additional redundant (i.e. correlated) features has the effect that the
# selected features vary depending on the cross-validation fold. The remaining
# features are non-informative as they are drawn at random.


from sklearn.datasets import make_classification

X, y = make_classification(
    n_samples=500,
    n_features=15,
    n_informative=3,
    n_redundant=2,
    n_repeated=0,
    n_classes=8,
    n_clusters_per_class=1,
    class_sep=0.8,
    random_state=0,
)

# Model training and selection
# ----------------------------
from sklearn.feature_selection import RFECV
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold

min_features_to_select = 1
clf = LogisticRegression(solver='liblinear', random_state=0)
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)

try:
    rfecv = RFECV(
        estimator=clf,
        step=1,
        cv=cv,
        scoring="accuracy",
        min_features_to_select=min_features_to_select,
        n_jobs=2,
    )
    rfecv.fit(X, y)
    print(f"Optimal number of features: {rfecv.n_features_}")
except Exception as e:
    print(f"An error occurred during fitting: {e}")

# Plot number of features VS. cross-validation scores
# ---------------------------------------------------
import matplotlib.pyplot as plt
import pandas as pd

try:
    cv_results = pd.DataFrame(rfecv.cv_results_)
    plt.figure()
    plt.xlabel("Number of features selected")
    plt.ylabel("Mean test accuracy")
    plt.errorbar(
        x=cv_results["n_features"],
        y=cv_results["mean_test_score"],
        yerr=cv_results["std_test_score"],
        fmt='-o',
        capsize=5,
    )
    plt.title("Recursive Feature Elimination \nwith correlated features")
    plt.show()
except Exception as e:
    print(f"An error occurred during plotting: {e}")
