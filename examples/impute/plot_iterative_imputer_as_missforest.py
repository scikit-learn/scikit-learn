"""
=============================================================
Replicating Functionality of missForest with IterativeImputer
=============================================================

There are many well-established imputation packages in the R data science
ecosystem: Amelia, mi, mice, missForest, and others.

missForest is popular, and turns out to be a particular instance of a class of
sequential imputation algorithms that can all be implemented with the
:class:`sklearn.impute.IterativeImputer` class, which is a strategy for
imputing missing values by modeling each feature with missing values as a
function of other features in a round-robin fashion. In the case of missForest,
the function is a Random Forest.

In this example we will demonstrate how to use
:class:`sklearn.impute.IterativeImputer` to replicate the functionality of
missForest.
"""
print(__doc__)

import numpy as np

from sklearn.datasets import load_diabetes
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import IterativeImputer

rng = np.random.RandomState(0)

# Load data
dataset = load_diabetes()
X_full, y_full = dataset.data, dataset.target
n_samples = X_full.shape[0]
n_features = X_full.shape[1]

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

# Random Forest predictor with default values according to missForest docs
predictor = RandomForestRegressor(n_estimators=100, max_features='sqrt')
imputer = IterativeImputer(n_iter=10, predictor=predictor)

# Impute missing values with IterativeImputer as missForest
X_imputed = imputer.fit_transform(X_missing)

# Compute RMSE of the imputed values
imp_missing_vals = X_imputed[np.where(missing_samples)[0], missing_features]
true_missing_vals = X_full[np.where(missing_samples)[0], missing_features]
rmse = np.sqrt(np.mean((true_missing_vals - imp_missing_vals)**2))
print('RMSE of IterativeImputer as missForest on the Diabetes Data:', rmse)
