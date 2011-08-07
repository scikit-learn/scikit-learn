"""
==========================
Cross validated Lasso path
==========================

Compute a 20-fold cross-validated :ref:`lasso` path to find the optimal value of
alpha.

This example uses two different algorithm to compute the Lasso path:
coordinate descent, as implemented by the LassoCV class, and Lars (least
angle regression) as implemented by the LassoLarsCV class. Both
algorithms give roughly the same results. They differ with regards to
their execution speed and sources of numerical errors.

Lars computes a path solution only for each kink in the path. As a
result, it is very efficient when there are only of few kinks, which is
the case if there are few features or samples. Also, it is able to
compute the full path without setting any meta parameter. On the
opposite, coordinate descent compute the path points on a pre-specified
grid (here we use the default). Thus it is more efficient if the number
of grid points is smaller than the number of kinks in the path. Such a
strategy can be interesting if the number of features is really large
and there are enough samples to select a large amount. In terms of
numerical errors, for heavily correlated variables, Lars will accumulate
more erros, while the coordinate descent algorithm will only sample the
path on a grid.

Note how the optimal value of alpha varies for each fold. This
illustrates why nested-cross validation is necessary when trying to
evaluate the performance of a method for which a parameter is chosen by
cross-validation: this choice of parameter may not be optimal for unseen
data.

Results obtained with LassoLarsIC which is based on AIC/BIC criteria
are also compared.
"""
print __doc__

# Author: Olivier Grisel, Gael Varoquaux, Alexandre Gramfort
# License: BSD Style.

import time

import numpy as np
import pylab as pl

from scikits.learn.linear_model import LassoCV, LassoLarsCV, LassoLarsIC
from scikits.learn import datasets

diabetes = datasets.load_diabetes()
X = diabetes.data
y = diabetes.target

rng = np.random.RandomState(42)
X = np.c_[X, rng.randn(X.shape[0], 14)]  # add some bad features

# normalize data as done by Lars to allow for comparison
X /= np.sqrt(np.sum(X ** 2, axis=0))

##############################################################################
# LassoLarsIC: least angle regression with BIC/AIC criterion

model_bic = LassoLarsIC(criterion='bic')
model_bic.fit(X, y)
alpha_bic_ = model_bic.alpha_

model_aic = LassoLarsIC(criterion='aic')
model_aic.fit(X, y)
alpha_aic_ = model_aic.alpha_

##############################################################################
# LassoCV: coordinate descent

# Compute paths
print "Computing regularization path using the coordinate descent lasso..."
t1 = time.time()
model = LassoCV(cv=20).fit(X, y)
t_lasso_cv = time.time() - t1

# Display results
m_log_alphas = -np.log10(model.alphas)

pl.figure()
ymin, ymax = 2300, 3800
pl.plot(m_log_alphas, model.mse_path_, ':')
pl.plot(m_log_alphas, model.mse_path_.mean(axis=-1), 'k',
        label='Average accross the folds', linewidth=2)
pl.axvline(-np.log10(model.alpha), linestyle='--', color='k',
           label='alpha CV')
pl.axvline(-np.log10(alpha_bic_), linestyle='--', color='r',
           label='alpha BIC')
pl.axvline(-np.log10(alpha_aic_), linestyle='--', color='g',
           label='alpha AIC')

pl.legend()

pl.xlabel('-log(lambda)')
pl.ylabel('Mean square error')
pl.title('Mean square error on each fold: coordinate descent (train time: %.2fs)' %
            t_lasso_cv)
pl.axis('tight')
pl.ylim(ymin, ymax)

##############################################################################
# LassoLarsCV: least angle regression

# Compute paths
print "Computing regularization path using the Lars lasso..."
t1 = time.time()
model = LassoLarsCV(cv=20).fit(X, y)
t_lasso_lars_cv = time.time() - t1

# Display results
m_log_alphas = -np.log10(model.cv_alphas)

pl.figure()
pl.plot(m_log_alphas, model.cv_mse_path_, ':')
pl.plot(m_log_alphas, model.cv_mse_path_.mean(axis=-1), 'k',
        label='Average accross the folds', linewidth=2)
pl.axvline(-np.log10(model.alpha), linestyle='--', color='k',
           label='alpha CV')
pl.axvline(-np.log10(alpha_bic_), linestyle='--', color='r',
           label='alpha BIC')
pl.axvline(-np.log10(alpha_aic_), linestyle='--', color='g',
           label='alpha AIC')
pl.legend()

pl.xlabel('-log(lambda)')
pl.ylabel('Mean square error')
pl.title('Mean square error on each fold: Lars (train time: %.2fs)' %
            t_lasso_lars_cv)
pl.axis('tight')
pl.ylim(ymin, ymax)

pl.show()
