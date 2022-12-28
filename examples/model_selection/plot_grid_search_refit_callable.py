"""
==================================================
Balance model complexity and cross-validated score
==================================================

This example balances model complexity and cross-validated score by
finding the most accurate model, while minimizing the number of PCA
components [1].

Based on the handwritten digits dataset and a Support Vector Classifier estimator,
the first figure shows the trade-off between cross-validated score and the number of
PCA components, when using a 1 standard-error rule to optimize model selection with
`GridSearchCV`. Under these constraints, the balanced case is when n_components=10 and
accuracy=0.88, which falls into the range within 1 standard error of the best accuracy
score.

Based on the Diabetes dataset and a Gradient Boosted Regressor, the second figure
shows the trade-off between cross-validated score and the number of tree estimators,
when using a 1 standard-error rule to optimize model selection with
`RandomizedSearchCV`. In this case, the balanced case is when n_estimators=100 and
R2=0.38, which falls into the range within 1 standard error of the best R2 score. This
example then closes with a brief exploratory analysis showing that the balanced case
lends to less test-set deviance when compared to that of the "unrazored",
highest-scoring case.

[1] Hastie, T., Tibshirani, R.,, Friedman, J. (2001). Model Assessment and
Selection. The Elements of Statistical Learning (pp. 219-260). New York,
NY, USA: Springer New York Inc..

"""

# Authors: Wenhao Zhang <wenhaoz@ucla.edu>, Derek Pisner <dpysalexander@gmail.com>

import numpy as np
import matplotlib.pyplot as plt

from sklearn import datasets
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV, Razors
from sklearn.pipeline import Pipeline
from sklearn.svm import LinearSVC
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split

# 1-SE rule - Digits dataset - GridSearchCV - PCA - LinearSVC - Classification
X, y = datasets.load_digits(return_X_y=True)
param_grid = {"reduce_dim__n_components": [6, 8, 10, 12, 14]}

pipe = Pipeline(
    [
        ("reduce_dim", PCA(random_state=42)),
        ("classify", LinearSVC(random_state=42, C=0.01)),
    ]
)

grid = GridSearchCV(
    pipe,
    cv=10,
    n_jobs=-1,
    param_grid=param_grid,
    scoring="accuracy",
    refit=Razors.simplify(param="reduce_dim__n_components", rule="se", sigma=1),
)

grid.fit(X, y)

n_components = grid.cv_results_["param_reduce_dim__n_components"]
test_scores = grid.cv_results_["mean_test_score"][
    ~grid.cv_results_["param_reduce_dim__n_components"].mask
]

plt.figure()
plt.bar(
    param_grid["reduce_dim__n_components"],
    grid.cv_results_["mean_test_score"],
    width=1.3,
    color="b",
)

lower = np.max(test_scores)
plt.axhline(np.max(test_scores), linestyle="--", color="y", label="Best score")
plt.axhline(lower, linestyle="--", color=".5", label="Best score - 1 SE")

plt.title("Balance model complexity and cross-validated score")
plt.xlabel("Number of PCA components used")
plt.ylabel("Digit classification accuracy")
plt.xticks(n_components.data.tolist())
plt.ylim((0.8, 1.0))
plt.legend(loc="upper left")

best_index_ = np.where(grid.cv_results_["mean_test_score"] == lower)[0][0]

print("The n_components selected is %d" % n_components.data[best_index_])
print(
    "The corresponding accuracy score is %.2f"
    % grid.cv_results_["mean_test_score"][best_index_]
)


# 1-SE rule - Diabetes dataset - RandomizedSearchCV - GradientBoostingRegressor -
# Regression
X, y = datasets.load_diabetes(return_X_y=True)

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.1, random_state=13
)

params = {
    "n_estimators": [25, 50, 75, 100, 150, 200, 250],
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
    cv=5,
    n_jobs=-1,
    scoring="r2",
    refit=Razors.simplify(param="n_estimators", rule="se", sigma=1),
)
grid.fit(X_train, y_train)

n_estimators = grid.cv_results_["param_n_estimators"]
test_scores = grid.cv_results_["mean_test_score"][
    ~grid.cv_results_["param_n_estimators"].mask
]

plt.bar(
    params["n_estimators"],
    grid.cv_results_["mean_test_score"],
    width=5,
    color="b",
)

lower = np.max(test_scores)
plt.axhline(lower, linestyle="--", color=".5", label="Best score - Wilcoxon 5%")
plt.axhline(
    np.max(grid.cv_results_["mean_test_score"]),
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

best_index_ = np.where(grid.cv_results_["mean_test_score"] == lower)[0][0]

print("The n_estimators selected is %d" % n_estimators.data[best_index_])
print(
    "The corresponding accuracy score is %.2f"
    % grid.cv_results_["mean_test_score"][best_index_]
)
plt.show()

# Show the effect on train-test deviance of using `Razors.simplify`
best_estimator_razor = grid.estimator
best_estimator_razor.n_estimators = params["n_estimators"][best_index_]
best_estimator_razor.fit(X_train, y_train)

# Compute test set deviance for razored and non-razored models
# razored model
test_score = np.zeros((params["n_estimators"][best_index_],), dtype=np.float64)
for i, y_pred in enumerate(best_estimator_razor.staged_predict(X_test)):
    test_score[i] = mean_squared_error(y_test, y_pred)

fig = plt.figure(figsize=(6, 6))
plt.subplot(1, 1, 1)
plt.title("Deviance of Razored Model")
plt.plot(
    np.arange(params["n_estimators"][best_index_]) + 1,
    best_estimator_razor.train_score_,
    "b-",
    label="Training Set Deviance",
)
plt.plot(
    np.arange(params["n_estimators"][best_index_]) + 1,
    test_score,
    "r-",
    label="Test Set Deviance",
)
plt.legend(loc="upper right")
plt.xlabel("Boosting Iterations")
plt.ylabel("Deviance")
fig.tight_layout()
plt.show()

best_estimator = grid.best_estimator_
best_estimator.fit(X_train, y_train)

# unrazored model
test_score = np.zeros((best_estimator.n_estimators,), dtype=np.float64)
for i, y_pred in enumerate(best_estimator.staged_predict(X_test)):
    test_score[i] = mean_squared_error(y_test, y_pred)

fig = plt.figure(figsize=(6, 6))
plt.subplot(1, 1, 1)
plt.title("Deviance of Non-Razored Model")
plt.plot(
    np.arange(best_estimator.n_estimators) + 1,
    best_estimator.train_score_,
    "b-",
    label="Training Set Deviance",
)
plt.plot(
    np.arange(best_estimator.n_estimators) + 1,
    test_score,
    "r-",
    label="Test Set Deviance",
)
plt.legend(loc="upper right")
plt.xlabel("Boosting Iterations")
plt.ylabel("Deviance")
fig.tight_layout()
plt.show()
