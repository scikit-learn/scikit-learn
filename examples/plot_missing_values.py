"""
======================================================
Imputing missing values before building an estimator
======================================================

This example shows that imputing the missing values can give better results
than discarding the samples containing any missing value. Imputing does not
always improve the predictions, so please check via cross-validation. Sometimes
dropping rows or using marker values is more effective.

Missing values can be replaced by the mean, the median, or the most frequent
value using the ``strategy`` hyper-parameter. The median is a more robust
estimator for data with high magnitude variables which could dominate results
(otherwise known as a 'long tail').

Another option is the MICE imputer. This uses round-robin linear regression,
treating every variable as an output in turn. The simple version implemented
assumes Gaussian output variables. If your output variables are obviously
non-Gaussian, consider transforming them to improve performance.

Script output::

    Results for the diabetes dataset:
    MSE with the entire dataset = 3354.15
    MSE without the samples containing missing values = 2968.98
    MSE after mean imputation of the missing values = 3507.77
    MSE after MICE imputation of the missing values = 3365.82

In this case, imputing helps the classifier get close to the original score
than using mean imputation.
"""

import numpy as np

from sklearn.datasets import load_diabetes
from sklearn.datasets import load_boston
from sklearn.ensemble import RandomForestRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import Imputer
from sklearn.preprocessing import MICEImputer
from sklearn.model_selection import cross_val_score

rng = np.random.RandomState(0)

def print_results(dataset):
    X_full, y_full = dataset.data, dataset.target
    n_samples = X_full.shape[0]
    n_features = X_full.shape[1]

    # Estimate the score on the entire dataset, with no missing values
    estimator = RandomForestRegressor(random_state=0, n_estimators=100)
    score = cross_val_score(estimator, X_full, y_full,
                            scoring='neg_mean_squared_error').mean() * -1
    print("MSE with the entire dataset = %.2f" % score)

    # Add missing values in 75% of the lines
    missing_rate = 0.75
    n_missing_samples = int(np.floor(n_samples * missing_rate))
    missing_samples = np.hstack((np.zeros(n_samples - n_missing_samples,
                                          dtype=np.bool),
                                 np.ones(n_missing_samples,
                                         dtype=np.bool)))
    rng.shuffle(missing_samples)
    missing_features = rng.randint(0, n_features, n_missing_samples)

    # Estimate the score without the lines containing missing values
    X_filtered = X_full[~missing_samples, :]
    y_filtered = y_full[~missing_samples]
    estimator = RandomForestRegressor(random_state=0, n_estimators=100)
    score = cross_val_score(estimator, X_filtered, y_filtered,
                            scoring='neg_mean_squared_error').mean() * -1
    print("MSE without the samples containing missing values = %.2f" % score)

    # Estimate the score after imputation (mean strategy) of the missing values
    X_missing = X_full.copy()
    X_missing[np.where(missing_samples)[0], missing_features] = 0
    y_missing = y_full.copy()
    estimator = Pipeline([("imputer", Imputer(missing_values=0,
                                              strategy="mean",
                                              axis=0)),
                          ("forest", RandomForestRegressor(random_state=0,
                                                           n_estimators=100))])
    score = cross_val_score(estimator, X_missing, y_missing,
                            scoring='neg_mean_squared_error').mean() * -1
    print("MSE after mean imputation of the missing values = %.2f" % score)

    # Estimate the score after imputation (MICE strategy) of the missing values
    estimator = Pipeline([("imputer", MICEImputer(missing_values=0,
                                                  random_state=0)),
                          ("forest", RandomForestRegressor(random_state=0,
                                                           n_estimators=100))])
    score = cross_val_score(estimator, X_missing, y_missing,
                            scoring='neg_mean_squared_error').mean() * -1
    print("MSE after MICE imputation of the missing values = %.2f" % score)


print("Results for the diabetes dataset:")
print_results(load_diabetes())

###########################################################################
# Note that MICE will not always be better than, e.g., simple mean imputation.
# To see an example of this, we swap in ``boston`` for ``diabetes``.
# Script output::
#
#    Results for the boston dataset:
#    MSE with the entire dataset = 28.95
#    MSE without the samples containing missing values = 31.93
#    MSE after mean imputation of the missing values = 29.41
#    MSE after MICE imputation of the missing values = 31.43

print("\nResults for the boston dataset:")
print_results(load_boston())
