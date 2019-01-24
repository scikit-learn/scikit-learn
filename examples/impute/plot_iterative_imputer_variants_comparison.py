"""
=========================================================
Imputing missing values with variants of IterativeImputer
=========================================================

The :class:`sklearn.impute.IterativeImputer` class is very flexible - it can be
used with a variety of predictors to do round-robin regression, treating every
variable as an output in turn.

In this example we compare some predictors for the purpose of missing feature
imputation with :class:`sklearn.imputeIterativeImputer`::

    :class:`sklearn.linear_model.BayesianRidge`: regularized linear regression
    :class:`sklearn.tree.DecisionTreeRegressor`: non-linear regression
    :class:`sklearn.neighbors.KNeighborsRegressor`: comparable to other KNN
    imputation approaches
    :class:`sklearn.ensemble.ExtraTreesRegressor`: similar to missForest in R

Of particular interest is the ability of
:class:`sklearn.impute.IterativeImputer` to mimic the behavior of missForest, a
popular imputation package for R. In this example, we have chosen to use
:class:`sklearn.ensemble.ExtraTreesRegressor` instead of
:class:`sklearn.ensemble.RandomForestRegressor` (as in missForest) due to its
increased speed.

Note that :class:`sklearn.neighbors.KNeighborsRegressor` is different from KNN
imputation, which learns from samples with missing values by using a distance
metric that accounts for missing values, rather than imputing them.

The goal is to compare different predictors to see which one is best for the
:class:`sklearn.impute.IterativeImputer` when using a
:class:`sklearn.linear_model.BayesianRidge` estimator on the California housing
dataset with a single value randomly removed from each row.

For this particular pattern of missing values we see that
:class:`sklearn.ensemble.ExtraTreesRegressor` and
:class:`sklearn.linear_model.BayesianRidge` give the best results.
"""
print(__doc__)

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from sklearn.datasets import fetch_california_housing
from sklearn.impute import SimpleImputer
from sklearn.impute import IterativeImputer
from sklearn.linear_model import BayesianRidge
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import cross_val_score

N_SPLITS = 5

rng = np.random.RandomState(0)

X_full, y_full = fetch_california_housing(return_X_y=True)
n_samples, n_features = X_full.shape

# Estimate the score on the entire dataset, with no missing values
br_estimator = BayesianRidge()
score_full_data = pd.DataFrame(
    cross_val_score(
        br_estimator, X_full, y_full, scoring='neg_mean_squared_error',
        cv=N_SPLITS
    )
)

# Add a single missing value to each row
X_missing = X_full.copy()
y_missing = y_full
missing_samples = np.arange(n_samples)
missing_features = rng.choice(n_features, n_samples, replace=True)
X_missing[missing_samples, missing_features] = np.nan

# Estimate the score after imputation (mean and median strategies) of the missing values
score_simple_imputer = pd.DataFrame()
for strategy in ('mean', 'median'):
    estimator = make_pipeline(
        SimpleImputer(missing_values=np.nan, strategy=strategy),
        br_estimator
    )
    score_simple_imputer[strategy] = cross_val_score(
        estimator, X_missing, y_missing, scoring='neg_mean_squared_error',
        cv=N_SPLITS
    )

# Estimate the score after iterative imputation of the missing values
# with different predictors
predictors = [
    BayesianRidge(),
    DecisionTreeRegressor(max_features='sqrt', random_state=0),
    KNeighborsRegressor(n_neighbors=15),
    ExtraTreesRegressor(n_estimators=10, n_jobs=-1, random_state=0)
]
score_iterative_imputer = pd.DataFrame()
for predictor in predictors:
    estimator = make_pipeline(
        IterativeImputer(random_state=0, predictor=predictor),
        br_estimator
    )
    score_iterative_imputer[predictor.__class__.__name__] = \
        cross_val_score(
            estimator, X_missing, y_missing, scoring='neg_mean_squared_error',
            cv=N_SPLITS
        )

scores = pd.concat(
    [score_full_data, score_simple_imputer, score_iterative_imputer],
    keys=['Original', 'SimpleImputer', 'IterativeImputer'], axis=1
)

labels = ['Full Data',
          'SimpleImputer w/ Mean Strategy',
          'SimpleImputer w/ Median Strategy',
          'IterativeImputer w/ BayesianRidge',
          'IterativeImputer w/ DecisionTreeRegressor',
          'IterativeImputer w/ KNeighborsRegressor',
          'IterativeImputer w/ ExtraTreesRegressor']

# plot boston results
fig, ax = plt.subplots(figsize=(13, 6))
means = -scores.mean()
errors = scores.std()
means.plot.barh(xerr=errors, ax=ax)
ax.set_title('California Housing Regression with Different Imputation Methods')
ax.set_xlabel('MSE (smaller is better)')
ax.set_yticks(np.arange(means.shape[0]))
ax.invert_yaxis()
ax.set_yticklabels(labels)
plt.tight_layout(pad=1)
plt.show()
