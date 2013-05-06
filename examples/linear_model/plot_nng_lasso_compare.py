#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
=========================================================
Nonnegative Garrote and Lars-path comparison
=========================================================
We replicate the comparison done by Ming Yuan and Yi Lin
in their paper
`On the Nonnegative Garrote Estimator
<www2.isye.gatech.edu/statistics/papers/05-25.pdf>`_
(page 11).

For different values of alpha over the generated data sets, we
count how many times the path of each contains at least one estimate
that covers the true model.
"""
print(__doc__)


# Code source: Jaques Grobler
# License: BSD

import numpy as np
import pylab as pl

from sklearn.linear_model import LinearRegression, Lasso, NonNegativeGarrote
from sklearn.linear_model import non_negative_garrote_path

from sklearn.utils import check_random_state
from sklearn.linear_model import lars_path

rng = check_random_state(None)

ng_path_correct = 0
lars_path_correct = 0
max_samples = 501

# True path for this experiment
coef = np.array([1, 1, 0])

# Lists for plotting the two techniques results
hits_lars = []
hits_ng = []

# For 4 different values of alpha
for alpha_val, fig_num in ((0.35, 1), (0.45, 2), (0.55, 3), (0.65, 4)):
    print 'for alpha = ', alpha_val
    # Set up plots
    pl.figure(fig_num, figsize=(5, 5))
    pl.clf
    pl.axis('tight')
    pl.title('alpha = %.2f' % alpha_val)
    pl.xlabel('Sample Size')
    pl.ylabel('Frequency of Selecting Correct Models')

    # Vary the sample size from 25 up until 500
    for sample_size in xrange(25, max_samples, 25):
        # Create 100 data sets to test
        for dataset_iter in xrange(0, 100):
            # Create a dataset
            # CHECK: Do example with 10 featues, three relevant, 3rd = X1 + X2
            X1 = rng.randn(sample_size)
            X2 = rng.randn(sample_size)
            X3 = (np.sqrt(1 - 2 * alpha_val ** 2) * rng.randn(sample_size)
                  + alpha_val * (X1 + X2))
            X = np.c_[X1, X2, X3]
            y = np.dot(X, coef) + rng.randn(sample_size)

            # Get the lasso path coefficients (using lars)
            _, _, lars_coefs = lars_path(X, y, method='lasso')

            # Get the non-negative garotte's coefficients
            ng_coefs, _ = non_negative_garrote_path(X, y, n_alphas=3)

            # Test if either model's solution path matches the original model
            # TODO: somehow pep8 these 2 long lines
            if np.any(np.all(ng_coefs.astype(np.bool) == coef.astype(np.bool)[:, np.newaxis], axis=0)):
                ng_path_correct = ng_path_correct + 1.0
            if np.any(np.all(lars_coefs.astype(np.bool) == coef.astype(np.bool)[:, np.newaxis], axis=0)):
                lars_path_correct = lars_path_correct + 1.0

            # Lasso_path version
            #if np.any(np.all(ng_coefs.astype(np.bool) == coef.astype(np.bool))):
            #        ng_path_correct = ng_path_correct + 1.0
            #if np.any(np.all(lars_coefs.astype(np.bool) == coef.astype(np.bool)[:, np.newaxis], axis=0)):
            #        lars_path_correct = lars_path_correct + 1.0

        hits_pers_lars = lars_path_correct / 100
        hits_lars.append(hits_pers_lars)
        lars_path_correct = 0
        hits_pers_ng = ng_path_correct / 100
        hits_ng.append(hits_pers_ng)
        ng_path_correct = 0

    pl.plot(xrange(25, max_samples, 25), hits_lars, 'r-')
    pl.plot(xrange(25, max_samples, 25), hits_ng, 'b-')
    pl.xlim([0, max_samples])
    pl.ylim([0, 1.1])
    hits_lars = []
    hits_ng = []

pl.show()
