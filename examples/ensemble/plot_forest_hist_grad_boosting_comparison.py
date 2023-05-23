"""
===============================================================
Comparing Random Forests and Histogram Gradient Boosting models
===============================================================

In this example we compare the performance of Random Forest (RF) and Histogram
Gradient Boosting (HGBT) models in terms of score and computation time for a
given dataset. Here we vary the parameters that control the number of trees
according to each estimator:

- `n_estimators` controls the number of trees in the forest. It's a fixed numer.
- `max_iter` is the the maximum number of iterations in a gradient boosting
  based model. The number of iterations corresponds to the number of trees for
  regression and binary classification problems. Furthermore, the actual number
  of trees required by the model depends on the stopping criteria.

HGBT use gradient boosting to iteratively improve the model's performance by
fitting each tree to the negative gradient of the loss function with respect to
the predicted value. RFs, on the other hand, are based on bagging and use a
majority vote to predict the outcome.

For more information on ensemble models, see the :ref:`User Guide <ensemble>`.
"""

# Author:  Arturo Amor <david-arturo.amor-quiroz@inria.fr>
# License: BSD 3 clause

# %%
# Load dataset
# ------------

from sklearn.datasets import fetch_california_housing

X, y = fetch_california_housing(return_X_y=True, as_frame=True)
n_samples, n_features = X.shape

# %%
# HGBT models may be faster than a parallelized RF in certain cases. For
# instance, HGBT can handle high-dimensional data more efficiently than RF due
# to its use of histogram-based binning, which reduces the number of splits that
# need to be evaluated when growing the trees. Also HGBT can be faster than RF
# for very large datasets, especially when the number of samples is larger than
# tens of thousands.
print(f"The dataset consists of {n_samples} samples and {n_features} features")

# %%
# Compute score and computation times
# -----------------------------------
#
# Notice that many parts of the implementation of
# :class:`~sklearn.ensemble.HistGradientBoostingClassifier` and
# :class:`~sklearn.ensemble.HistGradientBoostingRegressor` are parallelized by
# default.
#
# The implementation of :class:`~sklearn.ensemble.RandomForestRegressor` and
# :class:`~sklearn.ensemble.RandomForestClassifier` can also be run on multiple
# cores by using the `n_jobs` parameter, here set to 2 due to limitations on the
# documentation builder.

import joblib

N_CORES = joblib.cpu_count(only_physical_cores=True)
print(f"Number of physical cores: {N_CORES}")

# %%
# The other parameters of both models were tuned but the procedure is not shown
# here to keep the example simple.

from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GridSearchCV, KFold

models = {
    "Random Forest": RandomForestRegressor(
        min_samples_leaf=5, random_state=0, n_jobs=2
    ),
    "HistGradientBoosting": HistGradientBoostingRegressor(
        max_leaf_nodes=15, random_state=0
    ),
}
param_grids = {
    "Random Forest": {"n_estimators": [5, 10, 20, 50, 100]},
    "HistGradientBoosting": {"max_iter": [5, 10, 20, 50, 100, 200, 300]},
}
cv = KFold(n_splits=3, shuffle=True, random_state=0)
results = []

for name, model in models.items():
    grid_search = GridSearchCV(
        estimator=model,
        param_grid=param_grids[name],
        return_train_score=True,
        cv=cv,
    ).fit(X, y)
    result = {"model": name, "cv_results": grid_search.cv_results_}
    results.append(result)

# %%
# .. Note::
#  Tuning the `n_estimators` for RF generally results in a waste of computer
#  power. In practice one just needs to ensure that it is large enough so that
#  doubling its value does not lead to a significant improvement of the testing
#  score.
#
# Plot results
# ------------

import matplotlib.pyplot as plt

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))
for idx, result in enumerate(results):
    cv_results = result["cv_results"]

    param_name = list(param_grids[result["model"]].keys())[0]
    param_values = cv_results["param_" + param_name].compressed()
    train_scores = cv_results["mean_train_score"]
    train_score_err = cv_results["std_train_score"]
    test_scores = cv_results["mean_test_score"]
    test_score_err = cv_results["std_test_score"]
    train_times = cv_results["mean_fit_time"]
    train_time_err = cv_results["std_fit_time"]
    test_times = cv_results["mean_score_time"]
    test_time_err = cv_results["std_score_time"]

    axs[0, idx].errorbar(
        param_values, train_scores, yerr=train_score_err, label="Train score"
    )
    axs[0, idx].errorbar(
        param_values, test_scores, yerr=test_score_err, label="Test score"
    )
    axs[0, idx].set_xlabel(param_name)
    axs[0, idx].set_ylabel("Score")
    axs[0, idx].set_ylim(0.45, 0.95)
    axs[0, idx].set_title(result["model"])
    axs[0, idx].legend(loc="lower right")
    axs[1, idx].errorbar(
        param_values, train_times, yerr=train_time_err, label="Train time"
    )
    axs[1, idx].errorbar(
        param_values, test_times, yerr=test_time_err, label="Test time"
    )
    axs[1, idx].set_xlabel(param_name)
    axs[1, idx].set_ylabel("time (s)")
    axs[1, idx].set_yscale("log")
    axs[1, idx].set_ylim(0.001, 5)
    axs[1, idx].legend(loc="lower right")
plt.subplots_adjust(wspace=0.25)
plt.show()

# %%
# Both HGBT and RF models improve when increasing the number of trees in the
# ensemble. However, the scores reach a plateau where adding new trees just
# makes fitting and scoring slower.
#
# Unlike RF, HGBT models offer an early-stopping option (See
# :ref:`sphx_glr_auto_examples_ensemble_plot_gradient_boosting_early_stopping.py`)
# to avoid adding new unnecessary trees. Internally, the algorithm uses an
# out-of-sample set to compute the generalization performance of the model at
# each addition of a tree. Thus, if the generalization performance is not
# improving for more than `n_iter_no_change` iterations, it stops adding trees.
#
# Last but not least, in this example the training time of RF is much larger
# than the training time of HGBT, even for relatively low values of
# `n_estimators`. The reason is that boosting models rely on shallow trees,
# which predict faster. Nevertheless, the training and predicting time of RF can
# be reduced using the `n_jobs` parameter, as mentioned above.
#
# Overall, the performance of HGBT versus parallelized RF depends on the
# specific characteristics of the dataset and the modeling task. It's always a
# good idea to try both models and compare their performance on your specific
# problem to determine which model is the best fit.
