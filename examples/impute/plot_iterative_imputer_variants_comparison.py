"""
=========================================================
Imputing missing values with variants of IterativeImputer
=========================================================

The :class:`sklearn.impute.IterativeImputer` class is very flexible - it can be
used with a variety of predictors to do round-robin regression, treating every
variable as an output in turn.

In this example we compare some predictors for the purpose of missing feature
imputation with `IterativeImputer`::

    RidgeCV: regularized linear regression
    DecisionTreeRegressor: non-linear regression
    KNeighborsRegressor: comparable to other KNN imputation approaches
    ExtraTreesRegressor: similar to missForest in R

Of particular interest is the ability of ``IterativeImputer`` to mimic the
behavior of missForest, a popular imputation package for R. In this example,
we have chosen to use ``ExtraTreesRegressor`` instead of
``RandomForestRegressor`` (as in missForest) due to its increased speed.

The goal is to compare different predictors to see which one is best for
the `IterativeImputer` when using a ``BayesianRidge`` estimator on the
California housing dataset.

For this particular pattern of missing values we see that
``ExtraTreesRegressor`` gives the best results.
"""
print(__doc__)

import numpy as np
import matplotlib.pyplot as plt

from sklearn.datasets import fetch_california_housing
from sklearn.impute import SimpleImputer, IterativeImputer
from sklearn.linear_model import RidgeCV, BayesianRidge
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import cross_val_score

N_SPLITS = 5

rng = np.random.RandomState(0)

X_full, y_full = fetch_california_housing(return_X_y=True)
n_samples = X_full.shape[0]
n_features = X_full.shape[1]

# Estimate the score on the entire dataset, with no missing values
mses = np.zeros((7, N_SPLITS))
br_estimator = BayesianRidge()
mses[0, :] = cross_val_score(br_estimator, X_full, y_full,
                             scoring='neg_mean_squared_error',
                             cv=N_SPLITS)


# Add a single missing value to each row
X_missing = X_full.copy()
y_missing = y_full.copy()
n_missing_samples = int(np.floor(n_samples))
missing_samples = rng.choice(n_samples, n_missing_samples, replace=False)
missing_features = rng.choice(n_features, n_missing_samples, replace=True)
X_missing[missing_samples, missing_features] = np.nan

# Estimate the score after imputation (mean strategy) of the missing values
for i, strategy in enumerate(['mean', 'median']):
    estimator = make_pipeline(
        SimpleImputer(missing_values=np.nan, strategy=strategy),
        br_estimator
    )
    mses[i + 1, :] = cross_val_score(estimator, X_missing, y_missing,
                                     scoring='neg_mean_squared_error',
                                     cv=N_SPLITS)

# Estimate the score after iterative imputation of the missing values
# with different predictors
predictors = [
    RidgeCV(alphas=(1e-7, 0.01, 0.1, 1.0, 10.0)),
    DecisionTreeRegressor(random_state=0, max_features='sqrt'),
    KNeighborsRegressor(n_neighbors=15),
    ExtraTreesRegressor(n_estimators=10)
]

for i, predictor in enumerate(predictors):
    print(predictor)
    estimator = make_pipeline(
        IterativeImputer(random_state=0, predictor=predictor),
        br_estimator
    )
    mses[i + 3, :] = cross_val_score(estimator, X_missing, y_missing,
                                     scoring='neg_mean_squared_error',
                                     cv=N_SPLITS)

# Plot the results
x_labels = ['Full Data',
            'SimpleImputer w/ Mean Strategy',
            'SimpleImputer w/ Median Strategy',
            'IterativeImputer w/ RidgeCV',
            'IterativeImputer w/ DecisionTreeRegressor',
            'IterativeImputer w/ KNeighborsRegressor',
            'IterativeImputer w/ ExtraTreesRegressor']

# plot boston results
fig, ax = plt.subplots(figsize=(13, 6))
for i, j in enumerate(np.arange(mses.shape[0])):
    ax.barh(
        j,
        -np.mean(mses[j, :]),
        xerr=np.std(mses[j, :]),
        alpha=0.6,
        align='center'
    )

ax.set_title('California Housing Regression with Different Imputation Methods')
ax.set_xlabel('MSE (smaller is better)')
ax.set_yticks(np.arange(mses.shape[0]))
ax.invert_yaxis()
ax.set_yticklabels(x_labels)
plt.tight_layout(pad=1)
plt.show()
