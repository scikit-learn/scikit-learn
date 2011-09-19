"""
==================================================
Randomized Lasso: feature selection with Lasso
==================================================

"""
print __doc__

import pylab as pl

from sklearn.linear_model import RandomizedLasso
from sklearn.datasets import make_regression


################################################################################
# Generating simulated data with Gaussian weigthts
X, y, coef = make_regression(n_samples=40, n_features=100,
                             n_informative=10, coef=True, random_state=42)

################################################################################
# Fit the RandomizedLasso

################################################################################
# Plot the true weights, the estimated weights and the histogram of the
# weights
pl.figure(figsize=(6, 5))
pl.title("Weights of the model")
for index, alpha in enumerate((1, 5, 10, 15, 20, 'aic')):
    clf = RandomizedLasso(verbose=True, alpha=alpha, random_state=42)
    clf.fit(X, y)
    pl.plot(clf.scores_, ':', color=pl.cm.bone(index/6.),
                    label="alpha = %s" % alpha)
pl.plot(coef/coef.max(), 'g-', label="Ground truth coefficients")
pl.xlabel("Features")
#pl.ylabel("Values of the weights")
pl.legend(loc=1)

pl.show()
