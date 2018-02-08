"""
======================================================
Imputing missing values before building an estimator
======================================================

This example shows that imputing the missing values can give better
results than discarding the samples containing any missing value.
Imputing does not always improve the predictions, so please check via
cross-validation.  Sometimes dropping rows or using marker values is
more effective.

Missing values can be replaced by the mean, the median or the most frequent
value using the ``strategy`` hyper-parameter.
The median is a more robust estimator for data with high magnitude variables
which could dominate results (otherwise known as a 'long tail').

Script output::

  Score with the entire dataset = 0.88
  Score without the samples containing missing values = 0.68
  Score after imputation of the missing values (mean) = 0.88
  Score with NMF representation robust to missing values = 0.87

In this case, imputing helps the classifier get close to the original score.

"""
import numpy as np

from sklearn.datasets import fetch_olivetti_faces
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import Imputer
from sklearn.decomposition import NMF
from sklearn.model_selection import cross_val_score

print(__doc__)

rng = np.random.RandomState(0)

dataset = fetch_olivetti_faces()
X_full, y_full = dataset.data, dataset.target

n_samples = X_full.shape[0]
n_features = X_full.shape[1]

# Estimate the score on the entire dataset, with no missing values
estimator = RandomForestClassifier(random_state=0, n_estimators=100)
score = cross_val_score(estimator, X_full, y_full).mean()
print("Score with the entire dataset = %.2f" % score)

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
estimator = RandomForestClassifier(random_state=0, n_estimators=100)
score = cross_val_score(estimator, X_filtered, y_filtered).mean()
print("Score without the samples containing missing values = %.2f" % score)

# Estimate the score after imputation using a 'mean' strategy
X_missing = X_full.copy()
X_missing[np.where(missing_samples)[0], missing_features] = np.nan
y_missing = y_full.copy()
estimator = Pipeline([("imputer", Imputer(missing_values=np.nan,
                                          strategy="mean")),
                      ("forest", RandomForestClassifier(random_state=0,
                                                        n_estimators=100))])
score = cross_val_score(estimator, X_missing, y_missing).mean()
print("Score after imputation of the missing values (mean) = %.2f" % score)

# Estimate the score after imputation using non-negative matrix factorization
estimator = Pipeline([("nmf", NMF(solver='mu', init='random', n_components=19,
                                  random_state=0, max_iter=1000)),
                      ("forest", RandomForestClassifier(random_state=0,
                                                        n_estimators=100))])

score = cross_val_score(estimator, X_missing, y_missing).mean()
print("Score with NMF representation robust to missing values = %.2f" % score)
