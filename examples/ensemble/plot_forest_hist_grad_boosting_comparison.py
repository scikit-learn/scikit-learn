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
# documentation builder. See :ref:`parallelism` for more information.

import joblib

N_CORES = joblib.cpu_count(only_physical_cores=True)
print(f"Number of physical cores: {N_CORES}")

# %%
# The other parameters of both models were tuned but the procedure is not shown
# here to keep the example simple.

import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GridSearchCV, KFold

models = {
    "Random Forest": RandomForestRegressor(
        min_samples_leaf=5, random_state=0, n_jobs=2
    ),
    "HistGradientBoosting": HistGradientBoostingRegressor(
        max_leaf_nodes=15, random_state=0, early_stopping=False
    ),
}
param_grids = {
    "Random Forest": {"n_estimators": [10, 20, 50, 100]},
    "HistGradientBoosting": {"max_iter": [10, 20, 50, 100, 300, 500]},
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
    result = {"model": name, "cv_results": pd.DataFrame(grid_search.cv_results_)}
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
# We can use a `plotly.express.scatter
# <https://plotly.com/python-api-reference/generated/plotly.express.scatter.html>`_
# to visualize the trade-off between elapsed computing time and mean test score.
# Passing the cursor over a given point displays the corresponding parameters.
# Error bars correspond to one standard deviation as computed in the different
# folds of the cross-validation.

import plotly.express as px
import plotly.colors as colors
from plotly.subplots import make_subplots

fig = make_subplots(
    rows=1,
    cols=2,
    shared_yaxes=True,
    subplot_titles=["Train time vs score", "Predict time vs score"],
)
model_names = [result["model"] for result in results]
colors_list = colors.qualitative.Plotly * (
    len(model_names) // len(colors.qualitative.Plotly) + 1
)

for idx, result in enumerate(results):
    cv_results = result["cv_results"].round(3)
    model_name = result["model"]
    param_name = list(param_grids[model_name].keys())[0]
    cv_results[param_name] = cv_results["param_" + param_name]
    cv_results["model"] = model_name

    scatter_fig = px.scatter(
        cv_results,
        x="mean_fit_time",
        y="mean_test_score",
        error_x="std_fit_time",
        error_y="std_test_score",
        hover_data=param_name,
        color="model",
    )
    line_fig = px.line(
        cv_results,
        x="mean_fit_time",
        y="mean_test_score",
    )

    scatter_trace = scatter_fig["data"][0]
    line_trace = line_fig["data"][0]
    scatter_trace.update(marker=dict(color=colors_list[idx]))
    line_trace.update(line=dict(color=colors_list[idx]))
    fig.add_trace(scatter_trace, row=1, col=1)
    fig.add_trace(line_trace, row=1, col=1)

    scatter_fig = px.scatter(
        cv_results,
        x="mean_score_time",
        y="mean_test_score",
        error_x="std_score_time",
        error_y="std_test_score",
        hover_data=param_name,
    )
    line_fig = px.line(
        cv_results,
        x="mean_score_time",
        y="mean_test_score",
    )

    scatter_trace = scatter_fig["data"][0]
    line_trace = line_fig["data"][0]
    scatter_trace.update(marker=dict(color=colors_list[idx]))
    line_trace.update(line=dict(color=colors_list[idx]))
    fig.add_trace(scatter_trace, row=1, col=2)
    fig.add_trace(line_trace, row=1, col=2)

fig.update_layout(
    xaxis=dict(title="Train time (s) - lower is better"),
    yaxis=dict(title="Test R2 score - higher is better"),
    xaxis2=dict(title="Predict time (s) - lower is better"),
    legend=dict(x=0.72, y=0.05, traceorder="normal", borderwidth=1),
    title=dict(x=0.5, text="Speed-score trade-off of tree-based ensembles"),
)

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
# In this example HGBT performs uniformly better than RF in terms of test score
# and computing time (train and predict). Nevertheless, the training and
# predicting time of RF can be reduced using the `n_jobs` parameter, as
# mentioned above. Overall, the performance of HGBT versus parallelized RF
# depends on the specific characteristics of the dataset and the modeling task.
# It's always a good idea to try both models and compare their performance on
# your specific problem to determine which model is the best fit.
