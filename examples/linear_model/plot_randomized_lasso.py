"""
==============================================
Randomized Lasso: feature selection with Lasso
==============================================

"""
print __doc__

import pylab as pl
import numpy as np
from scipy.interpolate import interp1d

from sklearn.linear_model import RandomizedLasso, lars_path, LassoLarsCV
from sklearn.cross_validation import ShuffleSplit
from sklearn.datasets import load_diabetes
from sklearn.preprocessing import Scaler

###############################################################################
# Generating simulated data with Gaussian weigthts

diabetes = load_diabetes()
X = diabetes.data
y = diabetes.target

n_bad_features = 30
n_samples, n_good_features = X.shape

# Some noisy data not correlated
E = np.random.normal(size=(n_samples, n_bad_features))
# Add the noisy data to the informative features
X = np.hstack((X, E))
X = Scaler().fit_transform(X)

n_features = n_good_features + n_bad_features

###############################################################################
# Plot stability selection path

print "Computing regularization path using the LARS ..."

a = 0.2
n_resampling = 200
rng = np.random.RandomState(42)
xx_grid = np.linspace(0, 1, 100)
scores_path = np.zeros((n_features, 100))

for k in xrange(n_resampling):
    weights = 1. - a * rng.random_integers(0, 1, size=(n_features,))
    X_r = X * weights[np.newaxis, :]
    alphas, _, coefs = lars_path(X_r, y, method='lasso', verbose=False)

    xx = np.sum(np.abs(coefs.T), axis=1)
    xx /= xx[-1]

    interpolator = interp1d(xx, coefs)
    coefs_grid = interpolator(xx_grid)
    scores_path += (coefs_grid != 0.0)

scores_path /= n_resampling

pl.figure()
pl.plot(xx_grid, scores_path[:n_good_features].T, 'r')
pl.plot(xx_grid, scores_path[n_good_features:].T, 'k')
ymin, ymax = pl.ylim()
pl.xlabel('|coef| / max|coef|')
pl.ylabel('Proportion of times selected')
pl.title('Stability Scores Path')
pl.axis('tight')
pl.ylim([0, 1.05])

###############################################################################
# Plot the estimated stability scores for best cross-validated alpha

cv = ShuffleSplit(n_samples, n_iterations=50, test_fraction=0.2)
alpha_cv = LassoLarsCV(cv=cv).fit(X, y).alpha
clf = RandomizedLasso(verbose=True, alpha=alpha_cv, random_state=42, a=a)
clf.fit(X, y)

pl.figure()
pl.plot(clf.scores_)
pl.xlabel("Features")
pl.ylabel("Selection score")
pl.show()
