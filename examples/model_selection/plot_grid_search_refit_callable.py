"""
==================================================
Balance model complexity and cross-validated score
==================================================

This example balances model complexity and cross-validated score by
finding a decent accuracy within 1 standard deviation of the best accuracy
score while minimising the number of PCA components [1].

The figure shows the trade-off between cross-validated score and the number
of PCA components. The balanced case is when n_components=10 and accuracy=0.88,
which falls into the range within 1 standard deviation of the best accuracy
score.

[1] Hastie, T., Tibshirani, R.,, Friedman, J. (2001). Model Assessment and
Selection. The Elements of Statistical Learning (pp. 219-260). New York,
NY, USA: Springer New York Inc..

"""

# Author: Wenhao Zhang <wenhaoz@ucla.edu>

import numpy as np
import matplotlib.pyplot as plt

from sklearn.datasets import load_digits
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV, Razors
from sklearn.pipeline import Pipeline
from sklearn.svm import LinearSVC


def lower_bound(cv_results):
    """
    Calculate the lower bound within 1 standard deviation
    of the best `mean_test_scores`.

    Parameters
    ----------
    cv_results : dict of numpy(masked) ndarrays
        See attribute cv_results_ of `GridSearchCV`

    Returns
    -------
    float
        Lower bound within 1 standard deviation of the
        best `mean_test_score`.
    """
    best_score_idx = np.argmax(cv_results["mean_test_score"])

    return (
        cv_results["mean_test_score"][best_score_idx]
        - cv_results["std_test_score"][best_score_idx]
    )


pipe = Pipeline(
    [
        ("reduce_dim", PCA(random_state=42)),
        ("classify", LinearSVC(random_state=42, C=0.01)),
    ]
)

param_grid = {"reduce_dim__n_components": [6, 8, 10, 12, 14]}

grid = GridSearchCV(
    pipe,
    cv=10,
    n_jobs=-1,
    param_grid=param_grid,
    scoring="accuracy",
    refit=Razors.simplify(param="reduce_dim__n_components", rule="se", sigma=1),
)  # note that this would be similar to specifying rule="percentile", eta=0.68
X, y = load_digits(return_X_y=True)
grid.fit(X, y)

n_components = grid.cv_results_["param_reduce_dim__n_components"]
test_scores = grid.cv_results_["mean_test_score"]

plt.figure()
plt.bar(n_components, test_scores, width=1.3, color="b")

lower = lower_bound(grid.cv_results_)
plt.axhline(np.max(test_scores), linestyle="--", color="y", label="Best score")
plt.axhline(lower, linestyle="--", color=".5", label="Best score - 1 SE")

plt.title("Balance model complexity and cross-validated score")
plt.xlabel("Number of PCA components used")
plt.ylabel("Digit classification accuracy")
plt.xticks(n_components.data.tolist())
plt.ylim((0, 1.0))
plt.legend(loc="upper left")

best_index_ = grid.best_index_

print("The best_index_ is %d" % best_index_)
print("The n_components selected is %d" % n_components.data[best_index_])
print(
    "The corresponding accuracy score is %.2f"
    % grid.cv_results_["mean_test_score"][best_index_]
)
plt.show()


##
param_grid = {"reduce_dim__n_components": [6, 8, 10, 12, 14]}

grid = GridSearchCV(
    pipe,
    cv=10,
    n_jobs=-1,
    param_grid=param_grid,
    scoring="accuracy",
    refit=Razors.simplify(param="reduce_dim__n_components", rule="ranksum", alpha=0.05),
)
X, y = load_digits(return_X_y=True)
grid.fit(X, y)

n_components = grid.cv_results_["param_reduce_dim__n_components"]
test_scores = grid.cv_results_["mean_test_score"]

plt.figure()
plt.bar(n_components, test_scores, width=1.3, color="b")

lower = lower_bound(grid.cv_results_)
plt.axhline(np.max(test_scores), linestyle="--", color="y", label="Best score")
plt.axhline(lower, linestyle="--", color=".5", label="Best score - Wilcoxon 1%")

plt.title("Balance model complexity and cross-validated score")
plt.xlabel("Number of PCA components used")
plt.ylabel("Digit classification accuracy")
plt.xticks(n_components.data.tolist())
plt.ylim((0, 1.0))
plt.legend(loc="upper left")

best_index_ = grid.best_index_

print("The best_index_ is %d" % best_index_)
print("The n_components selected is %d" % n_components.data[best_index_])
print(
    "The corresponding accuracy score is %.2f"
    % grid.cv_results_["mean_test_score"][best_index_]
)
plt.show()
