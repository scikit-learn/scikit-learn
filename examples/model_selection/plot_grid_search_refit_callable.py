"""
==================================================
Balance model complexity and cross-validated score
==================================================

This example balances model complexity and cross-validated score by finding the most
accurate model, while minimizing the number of PCA components [1].

[1] Hastie, T., Tibshirani, R.,, Friedman, J. (2001). Model Assessment and
Selection. The Elements of Statistical Learning (pp. 219-260). New York,
NY, USA: Springer New York Inc..

"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np

from sklearn import datasets
from sklearn.decomposition import PCA
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import (
    FavorabilityRanker,
    GridSearchCV,
    RandomizedSearchCV,
    ScoreCutModelSelector,
    StandardErrorSlicer,
    WilcoxonSlicer,
    promote,
    train_test_split,
)
from sklearn.pipeline import Pipeline
from sklearn.svm import LinearSVC

"""
Based on the handwritten digits dataset and a Support Vector Classifier estimator,
the first figure shows the trade-off between cross-validated score and the number of
PCA components, when using a 1 standard-error rule to constrain model selection with
`GridSearchCV`. Under these constraints, the balanced case is when n_components=12 and
accuracy=0.90, which falls into the range within 1 standard error of the best accuracy
score.
"""

# 1-SE rule - Digits dataset - GridSearchCV - PCA - LinearSVC - Classification
X, y = datasets.load_digits(return_X_y=True)
param_grid = {"reduce_dim__n_components": [6, 8, 10, 12, 14]}

pipe = Pipeline(
    [
        ("reduce_dim", PCA(random_state=42)),
        ("classify", LinearSVC(random_state=42, C=0.01)),
    ]
)
favorability_rules = {
    "reduce_dim__n_components": (True, 1.0),  # Lower is simpler and more favorable
}
grid = GridSearchCV(
    pipe,
    cv=10,
    n_jobs=-1,
    param_grid=param_grid,
    scoring="accuracy",
    refit=promote(StandardErrorSlicer(sigma=1), FavorabilityRanker(favorability_rules)),
)

grid.fit(X, y)

n_components = grid.cv_results_["param_reduce_dim__n_components"]

original_best_index = np.argmax(grid.cv_results_["mean_test_score"])
original_best_score = np.max(grid.cv_results_["mean_test_score"])
simplest_index = grid.best_index_
simplest_score = grid.cv_results_["mean_test_score"][simplest_index]

plt.figure()
plt.bar(
    param_grid["reduce_dim__n_components"],
    grid.cv_results_["mean_test_score"],
    width=1.3,
    color="b",
)
plt.axhline(simplest_score, linestyle="--", color=".5", label="Best score - 1 SE")
plt.axhline(
    original_best_score,
    linestyle="--",
    color="y",
    label="Best score",
)

plt.title("Balance model complexity and cross-validated score")
plt.xlabel("Number of PCA components used")
plt.ylabel("Digit classification accuracy")
plt.xticks(n_components.data.tolist())
plt.ylim((0.8, 1.0))
plt.legend(loc="upper left")

plt.show()

"""
Based on the Diabetes dataset and a Gradient Boosted Regressor, the second figure
shows the trade-off between cross-validated score and the number of tree estimators,
when using a signed rank sum rule with an alpha-level of 0.05 to constrain model
selection with `RandomizedSearchCV`. Under these constraints, the balanced case is when
n_estimators=150 and R2=0.41, which is "not significantly different" from the best R2
score, but is nevertheless a more parsimonious model.

As an exploratory analysis that follows from the second example, we can observe that
the balanced case lends to less test-set deviance when compared to that of the more
complex highest-scoring case.
"""
# Signed Rank - Diabetes dataset - RandomizedSearchCV - GradientBoostingRegressor -
# Regression
plt.close("all")

X, y = datasets.load_diabetes(return_X_y=True)

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.1, random_state=13
)

params = {
    "n_estimators": [25, 50, 75, 100, 150, 200, 250],
}

favorability_rules = {
    "n_estimators": (True, 1.0),  # Lower is simpler and more favorable
}

reg = GradientBoostingRegressor(
    max_depth=4,
    min_samples_split=5,
    learning_rate=0.01,
    loss="squared_error",
    random_state=42,
)

grid = RandomizedSearchCV(
    reg,
    params,
    cv=10,
    n_jobs=-1,
    scoring="r2",
)
grid.fit(X_train, y_train)
n_estimators = grid.cv_results_["param_n_estimators"]
original_best_estimator = grid.best_estimator_
original_best_index = np.argmax(grid.cv_results_["mean_test_score"])
original_best_score = np.max(grid.cv_results_["mean_test_score"])
original_best_params = grid.cv_results_["params"][original_best_index]

plt.figure()
plt.bar(
    params["n_estimators"],
    grid.cv_results_["mean_test_score"],
    width=5,
    color="b",
)

plt.axhline(
    simplest_score, linestyle="--", color=".5", label="Best score - Wilcoxon 5%"
)
plt.axhline(
    original_best_score,
    linestyle="--",
    color="y",
    label="Best score",
)

plt.title("Balance model complexity and cross-validated score")
plt.xlabel("Number of estimators used")
plt.ylabel("Diabetes regression R2")
plt.xticks(
    np.arange(
        min(n_estimators.data.tolist()) - 25,
        max(n_estimators.data.tolist()) + 25,
        25.0,
    )
)
plt.ylim((0.3, 0.5))
plt.legend(loc="upper left")

# Compute test set deviance for original and promoted models
## Original model
original_best_estimator.fit(X_train, y_train)
original_test_score = np.zeros(
    (original_best_estimator.n_estimators,), dtype=np.float64
)
for i, y_pred in enumerate(original_best_estimator.staged_predict(X_test)):
    original_test_score[i] = mean_squared_error(y_test, y_pred)

fig1 = plt.figure(figsize=(6, 6))
plt.subplot(1, 1, 1)
plt.title("Deviance of Complex Model")
plt.plot(
    np.arange(original_best_estimator.n_estimators) + 1,
    original_best_estimator.train_score_,
    "b-",
    label="Training Set Deviance",
)
plt.plot(
    np.arange(original_best_estimator.n_estimators) + 1,
    original_test_score,
    "r-",
    label="Test Set Deviance",
)
plt.ylim((1500, 6000))
plt.yticks(np.arange(1500, 6000, 500))
plt.legend(loc="upper right")
plt.xlabel("Boosting Iterations")
plt.ylabel("Deviance")
fig1.tight_layout()

## Promoted model
ss = ScoreCutModelSelector(grid.cv_results_)
ss.fit(WilcoxonSlicer(alpha=0.05))
favorability_rules = {
    "reduce_dim__n_components": (True, 2.0),  # Lower is simpler and more favorable
    "classify__C": (False, 1.0),  # Lower is more complex and
    # less favorable
}
ss.transform(FavorabilityRanker(favorability_rules))
simplest_index = ss.favorable_best_index_
simplest_score = ss.favorable_best_score_
simplest_params = ss.favorable_best_params_
simplest_n_estimators = simplest_params["n_estimators"]

# Show the effect on train-test deviance of using constrained model selection
simplest_best_estimator = deepcopy(original_best_estimator)
simplest_best_estimator = original_best_estimator.set_params(**simplest_params)
simplest_best_estimator.fit(X_train, y_train)

simplest_test_score = np.zeros((simplest_n_estimators,), dtype=np.float64)
for i, y_pred in enumerate(simplest_best_estimator.staged_predict(X_test)):
    simplest_test_score[i] = mean_squared_error(y_test, y_pred)

fig2 = plt.figure(figsize=(6, 6))
plt.subplot(1, 1, 1)
plt.title("Deviance of Simpler Model")
plt.plot(
    np.arange(simplest_n_estimators) + 1,
    simplest_best_estimator.train_score_,
    "b-",
    label="Training Set Deviance",
)
plt.plot(
    np.arange(simplest_n_estimators) + 1,
    simplest_test_score,
    "r-",
    label="Test Set Deviance",
)
plt.ylim((1500, 6000))
plt.yticks(np.arange(1500, 6000, 500))
plt.legend(loc="upper right")
plt.xlabel("Boosting Iterations")
plt.ylabel("Deviance")
fig2.tight_layout()

plt.show()
