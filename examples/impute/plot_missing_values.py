"""
====================================================
Imputing missing values before building an estimator
====================================================

Missing values can be replaced by the mean, the median or the most frequent
value using the basic :class:`sklearn.impute.SimpleImputer`.

In this example we will investigate different imputation techniques on two
datasets: Diabetes dataset which is the set of parameteres collected from the
diabetes patients with aim to predict disease progression and California
Housing dataset for which the target is the median house value for California
districts.

Neither of those datasets has missing values. We will remove some of the values
and compare how will the results change if we use original data and the data
with imputed missing data by means of different techniques.

"""
print(__doc__)

# Authors: Maria Telenczuk  <https://github.com/maikia>
# License: BSD 3 clause

import numpy as np

rng = np.random.RandomState(0)
from sklearn.ensemble import RandomForestRegressor


N_SPLITS = 5
REGRESSOR = RandomForestRegressor(random_state=0)

###############################################################################
# Calculate the error
###############################################################################
#
# We are going to calculate the score for the imputers using negative mean
# square error. 

from sklearn.impute import MissingIndicator
from sklearn.model_selection import cross_val_score
from sklearn.pipeline import make_pipeline, make_union



def get_scores_for_imputer(imputer, X_missing, y_missing):
    estimator = make_pipeline(
        make_union(imputer, MissingIndicator(missing_values=0)),
        REGRESSOR)
    impute_scores = cross_val_score(estimator, X_missing, y_missing,
                                    scoring='neg_mean_squared_error',
                                    cv=N_SPLITS)
    return impute_scores

###############################################################################
#
###############################################################################
#The median is a more robust estimator for data with high magnitude variables
#which could dominate results (otherwise known as a 'long tail').
#
#With ``KNNImputer``, missing values can be imputed using the weighted
#or unweighted mean of the desired number of nearest neighbors.
#
# Another option is the :class:`sklearn.impute.IterativeImputer`. This uses
# round-robin linear regression, treating every variable as an output in
# turn. The version implemented assumes Gaussian (output) variables. If your
# features are obviously non-Normal, consider transforming them to look more
# Normal so as to potentially improve performance.
#
# In addition of using an imputing method, we can also keep an indication of the
# missing information using :func:`sklearn.impute.MissingIndicator` which might
#carry some information.
#

import matplotlib.pyplot as plt

# To use the experimental IterativeImputer, we need to explicitly ask for it:
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.datasets import load_diabetes


from sklearn.impute import (
    SimpleImputer, KNNImputer, IterativeImputer)


def get_results(dataset):
    X_full, y_full = dataset.data, dataset.target

    # As California dataset is quite large [(20640, 8)] to speed up the
    # calculations we will limit number of entries of both datasets to 440,
    # however feel free to use the whole datasets

    X_full = X_full[:440]
    y_full = y_full[:440]

    n_samples = X_full.shape[0]
    n_features = X_full.shape[1]

    # Estimate the score on the entire dataset, with no missing values
    full_scores = cross_val_score(REGRESSOR, X_full, y_full,
                                  scoring='neg_mean_squared_error',
                                  cv=N_SPLITS)

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
    X_missing[np.where(missing_samples)[0], missing_features] = 0
    y_missing = y_full.copy()

    # Estimate the score after replacing missing values by 0
    imputer = SimpleImputer(missing_values=0,
                            strategy='constant',
                            fill_value=0)
    zero_impute_scores = get_scores_for_imputer(imputer, X_missing, y_missing)

    # Estimate the score after imputation (mean strategy) of the missing values
    imputer = SimpleImputer(missing_values=0, strategy="mean")
    mean_impute_scores = get_scores_for_imputer(imputer, X_missing, y_missing)

    # Estimate the score after kNN-imputation of the missing values
    imputer = KNNImputer(missing_values=0)
    knn_impute_scores = get_scores_for_imputer(imputer, X_missing, y_missing)

    # Estimate the score after iterative imputation of the missing values
    imputer = IterativeImputer(missing_values=0,
                               random_state=0,
                               n_nearest_features=5,
                               sample_posterior=True)
    iterative_impute_scores = get_scores_for_imputer(imputer,
                                                     X_missing,
                                                     y_missing)

    return ((full_scores.mean(), full_scores.std()),
            (zero_impute_scores.mean(), zero_impute_scores.std()),
            (mean_impute_scores.mean(), mean_impute_scores.std()),
            (knn_impute_scores.mean(), knn_impute_scores.std()),
            (iterative_impute_scores.mean(), iterative_impute_scores.std()))


###############################################################################
# Download the data
###############################################################################
#
from sklearn.datasets import fetch_california_housing

results_diabetes = np.array(get_results(load_diabetes()))
mses_diabetes = results_diabetes[:, 0] * -1
stds_diabetes = results_diabetes[:, 1]

results_california = np.array(get_results(fetch_california_housing()))
mses_california = results_california[:, 0] * -1
stds_california = results_california[:, 1]

n_bars = len(mses_diabetes)
xval = np.arange(n_bars)

x_labels = ['Full data',
            'Zero imputation',
            'Mean Imputation',
            'KNN Imputation',
            'Iterative Imputation']
colors = ['r', 'g', 'b', 'orange', 'black']

# plot diabetes results
plt.figure(figsize=(12, 6))
ax1 = plt.subplot(121)
for j in xval:
    ax1.barh(j, mses_diabetes[j], xerr=stds_diabetes[j],
             color=colors[j], alpha=0.6, align='center')

ax1.set_title('Imputation Techniques with Diabetes Data')
ax1.set_xlim(left=np.min(mses_diabetes) * 0.9,
             right=np.max(mses_diabetes) * 1.1)
ax1.set_yticks(xval)
ax1.set_xlabel('MSE')
ax1.invert_yaxis()
ax1.set_yticklabels(x_labels)

# plot california results
ax2 = plt.subplot(122)
for j in xval:
    ax2.barh(j, mses_california[j], xerr=stds_california[j],
             color=colors[j], alpha=0.6, align='center')

# plot Ames results

ax2.set_title('Imputation Techniques with California Data')
ax2.set_yticks(xval)
ax2.set_xlabel('MSE')
ax2.invert_yaxis()
ax2.set_yticklabels([''] * n_bars)

plt.show()
