"""
==============================================
Randomized Lasso: feature selection with Lasso
==============================================

Performs feature scoring and selection using a randomized linear
model (Lasso)

"""
print __doc__

import pylab as pl
import numpy as np

from sklearn.linear_model import RandomizedLasso, lasso_stability_path, \
                                 LassoLarsCV
from sklearn.cross_validation import ShuffleSplit
from sklearn.datasets import load_diabetes
from sklearn.feature_selection import f_regression
from sklearn.preprocessing import Scaler

###############################################################################
# Load diabetes data and add noisy features
diabetes = load_diabetes()
X = diabetes.data
y = diabetes.target
n_bad_features = 15
n_samples, n_good_features = X.shape
n_features = n_good_features + n_bad_features
# Noisy data to the informative features
X = np.hstack((X, np.random.normal(size=(n_samples, n_bad_features))))
X = Scaler().fit_transform(X)

###############################################################################
# Plot stability selection path

print "Computing stability path using the LARS ..."

scaling = 0.3
coef_grid, scores_path = lasso_stability_path(X, y, scaling=scaling,
                                      random_state=42, sample_fraction=0.75)

pl.figure()
hg = pl.plot(coef_grid, scores_path[:n_good_features].T, 'r')
hb = pl.plot(coef_grid, scores_path[n_good_features:].T, 'k')
ymin, ymax = pl.ylim()
pl.xlabel('|coef| / max|coef|')
pl.ylabel('Proportion of times selected')
pl.title('Stability Scores Path')
pl.axis('tight')
pl.legend((hg[0], hb[0]), ('good features', 'bad features'),
          loc='lower right')
pl.ylim([0, 1.05])

###############################################################################
# Plot the estimated stability scores for best cross-validated alpha

cv = ShuffleSplit(n_samples, n_iterations=50, test_fraction=0.2)
alpha_cv = LassoLarsCV(cv=cv).fit(X, y).alpha
clf = RandomizedLasso(verbose=False, alpha=alpha_cv, random_state=42,
                      scaling=scaling)
clf.fit(X, y)

# compare with F-score
F, _ = f_regression(X, y)

pl.figure()
pl.plot(clf.scores_, label="Stabiliy selection score")
pl.plot(F / np.max(F), label="(scaled) F-score")
pl.xlabel("Features")
pl.ylabel("Score")
pl.legend()
pl.show()
