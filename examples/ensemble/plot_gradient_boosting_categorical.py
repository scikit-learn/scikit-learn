"""
========================================
Categorical Support in Gradient Boosting
========================================

.. currentmodule:: sklearn

In this example, we will compare the performance of
:class:`~ensemble.HistGradientBoostingRegressor` using one hot encoding
and with native categorical support.

We will work with the Ames Lowa Housing dataset which consists of numerical
and categorical features, where the houses' sales prices is the target.
"""
##############################################################################
# Load Ames Housing dataset
# -------------------------
# First, we load the ames housing data as a pandas dataframe. The features
# are either categorical or numerical:
print(__doc__)

from sklearn.datasets import fetch_openml

X, y = fetch_openml(data_id=41211, as_frame=True, return_X_y=True)

n_features = X.shape[1]
n_categorical_features = (X.dtypes == 'category').sum()
n_numerical_features = (X.dtypes == 'float').sum()
print(f"Number of features: {X.shape[1]}")
print(f"Number of categorical featuers: {n_categorical_features}")
print(f"Number of numerical featuers: {n_numerical_features}")

##############################################################################
# Create gradient boosting estimator with feature dropped
# ------------------------------------------------------------------
# As a baseline, we create a estimator where the categorical features are
# dropped:

from sklearn.experimental import enable_hist_gradient_boosting  # noqa
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.pipeline import make_pipeline
from sklearn.compose import make_column_transformer
from sklearn.compose import make_column_selector

preprocessor = make_column_transformer(
    ('drop', make_column_selector(dtype_include='category')),
    remainder='passthrough')
hist_dropped = make_pipeline(
    preprocessor, HistGradientBoostingRegressor(random_state=42))

##############################################################################
# Create gradient boosting estimator with one hot encoding
# --------------------------------------------------------
# Next, we create a pipeline that will one hot encode the categorical features
# and let rest of the numerical data to passthrough:

from sklearn.preprocessing import OneHotEncoder

preprocessor = make_column_transformer(
    (OneHotEncoder(sparse=False, handle_unknown='ignore'),
     make_column_selector(dtype_include='category')),
    remainder='passthrough')

hist_one_hot = make_pipeline(preprocessor,
                             HistGradientBoostingRegressor(random_state=42))

##############################################################################
# Create gradient boosting estimator with native categorical support
# ------------------------------------------------------------------
# The :class:`~ensemble.HistGradientBoostingRegressor` has native support
# for categorical features using the `categorical_features` parameter:

from sklearn.preprocessing import OrdinalEncoder

preprocessor = make_column_transformer(
    (OrdinalEncoder(handle_unknown='use_encoded_value', unknown_value=-1),
     make_column_selector(dtype_include='category')),
    remainder='passthrough')
categorical_mask = ([True] * n_categorical_features +
                    [False] * n_numerical_features)
hist_native = make_pipeline(
    preprocessor, HistGradientBoostingRegressor(
        random_state=42, categorical_features=categorical_mask))

##############################################################################
# Train the models with cross-validation
# --------------------------------
# Finally, we train the models using cross validation. Here we compare the
# models performance in terms of :func:`~metrics.r2_score` and fit times. We
# show that fit times are faster with native categorical support and that the
# test scores and scores times are comparable:

from sklearn.model_selection import cross_validate
import matplotlib.pyplot as plt
import numpy as np

dropped_result = cross_validate(hist_dropped, X, y, cv=3)
native_result = cross_validate(hist_native, X, y, cv=3)
one_hot_result = cross_validate(hist_one_hot, X, y, cv=3)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 8))

plot_info = [('fit_time', 'Fit times (s)', ax1, None),
             ('test_score', 'Test Scores (r2 score)', ax2, (0.5, 1.0))]

x, width = np.arange(3), 0.9
for key, title, ax, y_limit in plot_info:
    items = [native_result[key], dropped_result[key], one_hot_result[key]]
    ax.bar(x, [np.mean(item) for item in items],
           width, yerr=[np.std(item) for item in items],
           color=['b', 'r', 'g'])
    ax.set(xlabel='Type of model', title=title, xticks=x,
           xticklabels=['Native', "Dropped", "One Hot"], ylim=y_limit)
plt.show()
