"""
====================================================
Imputing missing values before building an estimator
====================================================

Missing values can be replaced by the mean, the median or the most frequent
value using the basic ``Imputer``.
The median is a more robust estimator for data with high magnitude variables
which could dominate results (otherwise known as a 'long tail').

Another option is the MICE imputer. This uses round-robin linear regression,
treating every variable as an output in turn. The version implemented assumes
Gaussian (output) variables. If your features are obviously non-Normal,
consider transforming them to look more Normal so as to improve performance.
"""

import numpy as np
import matplotlib.pyplot as plt

from sklearn.datasets import load_diabetes
from sklearn.datasets import load_boston
from sklearn.ensemble import RandomForestRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import Imputer
from sklearn.preprocessing import MICEImputer
from sklearn.model_selection import cross_val_score

rng = np.random.RandomState(0)


def get_results(dataset):
    X_full, y_full = dataset.data, dataset.target
    n_samples = X_full.shape[0]
    n_features = X_full.shape[1]

    # Estimate the score on the entire dataset, with no missing values
    estimator = RandomForestRegressor(random_state=0, n_estimators=100)
    full_score = cross_val_score(estimator, X_full, y_full,
                                 scoring='neg_mean_squared_error').mean()

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
    subset_score = cross_val_score(estimator, X_filtered, y_filtered,
                                   scoring='neg_mean_squared_error').mean()

    # Estimate the score after imputation (mean strategy) of the missing values
    X_missing = X_full.copy()
    X_missing[np.where(missing_samples)[0], missing_features] = 0
    y_missing = y_full.copy()
    estimator = Pipeline([("imputer", Imputer(missing_values=0,
                                              strategy="mean",
                                              axis=0)),
                          ("forest", RandomForestRegressor(random_state=0,
                                                           n_estimators=100))])
    mean_impute_score = cross_val_score(estimator, X_missing, y_missing,
                                        scoring='neg_mean_squared_error'
                                        ).mean()

    # Estimate the score after imputation (MICE strategy) of the missing values
    estimator = Pipeline([("imputer", MICEImputer(missing_values=0,
                                                  random_state=0)),
                          ("forest", RandomForestRegressor(random_state=0,
                                                           n_estimators=100))])
    mice_impute_score = cross_val_score(estimator, X_missing, y_missing,
                                        scoring='neg_mean_squared_error'
                                        ).mean()

    return full_score, subset_score, mean_impute_score, mice_impute_score


mses_diabetes = np.array(get_results(load_diabetes())) * -1
mses_boston = np.array(get_results(load_boston())) * -1

plt.figure(figsize=(12, 6))
xval = np.arange(len(mses_diabetes))
labels = {0: ('MSE with entire dataset', 'r'),
          1: ('MSE excluding samples containing missing values', 'g'),
          2: ('MSE after mean imputation of the missing values', 'b'),
          3: ('MSE after MICE imputation of the missing values', 'orange')}

# plot diabetes results
ax1 = plt.subplot(121)
for j in range(len(xval)):
    label, color = labels[xval[j]]
    ax1.bar(xval[j], mses_diabetes[j], width=0.8, color=color, alpha=0.6,
            align='center', label=label)
ax1.legend(loc='upper left')
ax1.set_ylabel('Mean Squared Error')
ax1.set_title('Feature Selection Techniques with Diabetes Data')
ax1.set_ylim(bottom=np.min(mses_diabetes) * 0.9,
             top=np.max(mses_diabetes) * 1.15)

# plot boston results
ax2 = plt.subplot(122)
for j in range(len(xval)):
    label, color = labels[xval[j]]
    ax2.bar(xval[j], mses_boston[j], width=0.8, color=color, alpha=0.6,
            align='center', label=label)
ax2.legend(loc='upper left')
ax2.set_ylabel('Mean Squared Error')
ax2.set_title('Feature Selection Techniques with Boston Data')
ax2.set_ylim(bottom=np.min(mses_boston) * 0.9,
             top=np.max(mses_boston) * 1.15)

plt.show()
