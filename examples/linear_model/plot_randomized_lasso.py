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

from sklearn import datasets
from sklearn.feature_selection import f_regression
from sklearn.preprocessing import Scaler

###############################################################################
# Load diabetes data and add noisy features
n_samples = 40
n_features = 100
n_good_features = 5
X, y, coef = datasets.make_regression(n_samples=n_samples,
                    n_features=n_features, n_informative=n_good_features,
                    coef=True, random_state=42)
X = Scaler().fit_transform(X)

###############################################################################
# Plot stability selection path
scaling = .3
alpha_grid, scores_path = lasso_stability_path(X, y, scaling=scaling,
                                      random_state=42)

pl.figure()
hg = pl.plot((alpha_grid[1:])**(.333), scores_path[coef != 0].T[1:], 'r')
hb = pl.plot((alpha_grid[1:])**(.333), scores_path[coef == 0].T[1:], 'k')
ymin, ymax = pl.ylim()
pl.xlabel(r'$(\alpha / \alpha_{max})^{1/3}$')
pl.ylabel('Proportion of times selected')
pl.title('Stability Scores Path')
pl.axis('tight')
pl.legend((hg[0], hb[0]), ('good features', 'bad features'),
          loc='lower right')
pl.ylim([0, 1.05])
#pl.xlim(0, 5)

###############################################################################
# Plot the estimated stability scores for best cross-validated alpha

# First find the best alpha:
cv = ShuffleSplit(n_samples, n_iterations=10, test_fraction=0.25)
ms = LassoLarsCV().fit(X, y)
alpha_cv = ms.alpha

# Then run the RandomizedLasso
clf = RandomizedLasso(verbose=False, alpha=alpha_cv, random_state=42,
                      scaling=scaling)
clf.fit(X, y)

F, _ = f_regression(X, y)  # compare with F-score

pl.figure()
score_plot = pl.plot(clf.scores_ / np.max(clf.scores_))
f_plot = pl.plot(F / np.max(F))
ground_truth = pl.plot(np.where(coef != 0), .5, 'ro')
pl.xlabel("Features")
pl.ylabel("Score")
pl.legend((score_plot, f_plot, ground_truth[0]),
            ("Stability-selection score", "(scaled) F-score",
             "Ground truth"))
pl.show()
