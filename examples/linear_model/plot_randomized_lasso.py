"""
==================================================
Randomized Lasso: feature selection with Lasso
==================================================

"""
print __doc__

import pylab as pl
import numpy as np

from sklearn.linear_model import RandomizedLasso, LassoLarsIC, LassoLarsCV
from sklearn.datasets import make_regression


################################################################################
# Generating simulated data with Gaussian weigthts
X, y, coef = make_regression(n_samples=40, n_features=100,
                             n_informative=10, coef=True, random_state=42)

################################################################################
# Plot the true weights, the estimated weights and the histogram of the
# weights
pl.figure(figsize=(6, 5))
pl.title("Weights of the model")
scores = list()
scores.append(coef/coef.max())
n_alphas = 12
alphas = np.logspace(-1, 1.5, n_alphas)[::-1]
for alpha in alphas:
    clf = RandomizedLasso(verbose=True, alpha=alpha, random_state=42)
    clf.fit(X, y)
    scores.append(clf.scores_)
scores = np.array(scores)
pl.imshow(np.rot90(scores), cmap=pl.cm.spectral,
          interpolation='nearest', aspect='auto')
#pl.plot(coef/coef.max(), 'g-', label="Ground truth coefficients")
alpha_aic = LassoLarsIC(criterion='aic').fit(X, y).alpha_
alpha_bic = LassoLarsIC(criterion='bic').fit(X, y).alpha_
alpha_cv = LassoLarsCV().fit(X, y).alpha
pl.xlabel("alpha (AIC=%.2e, BIC=%.2e, CV=%.2e)" %
          (alpha_aic, alpha_bic, alpha_cv))
pl.ylabel("Features")

pl.xticks(np.arange(0, n_alphas+1),
          ['Truth', ] + ['%.1e' % a for a in alphas])
pl.colorbar()
#pl.legend(loc=1)

pl.show()

