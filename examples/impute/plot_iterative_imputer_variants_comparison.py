"""
=========================================================
Imputing missing values with variants of IterativeImputer
=========================================================

The :class:`sklearn.impute.IterativeImputer` class is very flexible - it can be
used with a variety of predictors to do round-robin regression, treating every
variable as an output in turn.

In this example we compare some predictors for the purpose of missing feature
imputation with `IterativeImputer`::

    RidgeCV: default
    HuberRegressor: robust linear regression to reduce the impact of outliers
    DecisionTreeRegressor: non-linear regression
    RandomForestRegressor: equivalent to missForest in R
    KNeighborsRegressor: comparable to other KNN imputation approaches

The goal is to compare different predictors to see which one is best for
the `IterativeImputer` when using a ``RandomForestRegressor`` estimator on the
Boston dataset.

For the Boston dataset, the ``HuberRegressor`` produces results that are on
average superior to even having the full dataset. We also see that using other
predictors results in an imputer that is worse than using ``SimpleImputer``
with the ``mean`` strategy.
"""
print(__doc__)

import numpy as np
import matplotlib.pyplot as plt

from sklearn.datasets import load_boston
from sklearn.impute import SimpleImputer, IterativeImputer
from sklearn.linear_model import RidgeCV, HuberRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import cross_val_score

rng = np.random.RandomState(0)

X_full, y_full = load_boston(return_X_y=True)
n_samples = X_full.shape[0]
n_features = X_full.shape[1]

# Estimate the score on the entire dataset, with no missing values
rf_estimator = RandomForestRegressor(random_state=0, n_estimators=100)
full_scores = cross_val_score(rf_estimator, X_full, y_full,
                              scoring='neg_mean_squared_error',
                              cv=5)
mses_boston = [-full_scores.mean()]
stds_boston = [full_scores.std()]

# Add missing values in 75% of the lines
missing_rate = 0.75
n_missing_samples = int(np.floor(n_samples * missing_rate))
missing_samples = np.hstack((np.zeros(n_samples - n_missing_samples,
                                      dtype=np.bool),
                             np.ones(n_missing_samples,
                                     dtype=np.bool)))
rng.shuffle(missing_samples)
missing_features = rng.randint(0, n_features, n_missing_samples)
X_missing = X_full.copy()
X_missing[np.where(missing_samples)[0], missing_features] = np.nan
y_missing = y_full.copy()

# Estimate the score after imputation (mean strategy) of the missing values
for strategy in ['mean', 'median']:
    estimator = make_pipeline(
        SimpleImputer(missing_values=np.nan, strategy=strategy),
        rf_estimator
    )
    mean_impute_scores = cross_val_score(estimator, X_missing, y_missing,
                                         scoring='neg_mean_squared_error',
                                         cv=5)
    mses_boston.append(-mean_impute_scores.mean())
    stds_boston.append(mean_impute_scores.std())

# Estimate the score after iterative imputation of the missing values
# with different predictors
predictors = [
    RidgeCV(alphas=(1e-7, 0.01, 0.1, 1.0, 10.0)),
    HuberRegressor(),
    DecisionTreeRegressor(random_state=0, max_features='sqrt'),
    RandomForestRegressor(random_state=0,
                          n_estimators=100,
                          max_features='sqrt'),
    KNeighborsRegressor(n_neighbors=15)
]

for predictor in predictors:
    estimator = make_pipeline(
        IterativeImputer(random_state=0, predictor=predictor),
        rf_estimator
    )
    pred_scores = cross_val_score(estimator, X_missing, y_missing,
                                  scoring='neg_mean_squared_error',
                                  cv=5)
    mses_boston.append(-pred_scores.mean())
    stds_boston.append(pred_scores.std())

# Plot the results
x_labels = ['Full Data',
            'SimpleImputer w/ Mean Strategy',
            'SimpleImputer w/ Median Strategy',
            'IterativeImputer w/ RidgeCV',
            'IterativeImputer w/ HuberRegressor',
            'IterativeImputer w/ DecisionTreeRegressor',
            'IterativeImputer w/ RandomForestRegressor',
            'IterativeImputer w/ KNeighborsRegressor']

# plot boston results
fig, ax = plt.subplots(figsize=(10, 6))
for i, j in enumerate(np.arange(len(mses_boston))):
    ax.barh(j, mses_boston[j], xerr=stds_boston[j], alpha=0.6, align='center')

ax.set_title('Boston Data Regression MSE With Different Imputation Methods')
ax.set_xlabel('MSE')
ax.set_yticks(np.arange(len(mses_boston)))
ax.invert_yaxis()
ax.set_yticklabels(x_labels)
plt.show()
