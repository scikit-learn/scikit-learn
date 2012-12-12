# http://www2.isye.gatech.edu/statistics/papers/05-25.pdf , page 11
#
#imports
import numpy as np
import pylab as pl

from sklearn.linear_model import LinearRegression, Lasso, NonNegativeGarrote
from sklearn.linear_model import non_negative_garotte_path

from sklearn.utils import check_random_state
from sklearn.linear_model import lars_path

print 'running nngarrote.py - __main__'
rng = check_random_state(None)

ng_path_correct = 0
lars_path_correct = 0
max_samples = 501

# true path
coef = np.array([1, 1, 0])

# lists for plotting the two techniques results
hits_lars = []
hits_ng = []

# for 4 different values of alpha
for alpha_val, fig_num in ((0.35, 1), (0.45, 2), (0.55, 3), (0.65, 4)):
    print 'for alpha = ', alpha_val
    # set up plots
    pl.figure(fig_num, figsize=(5, 5))
    pl.clf
    pl.axis('tight')
    pl.title('alpha = %.2f' % alpha_val )
    pl.xlabel('Sample Size')
    pl.ylabel('Frequency of Selecting Correct Models')

    # vary the sample size from 25 up until 500
    for sample_size in xrange(25, max_samples, 25):
        # create 100 data sets to test
        for dataset_iter in xrange(0, 100):
            # create a dataset
            # CHECK: Do example with 10 featues, three relevant, 3rd = X1 + X2
            X1 = rng.randn(sample_size)
            X2 = rng.randn(sample_size)
            X3 = (np.sqrt(1 - 2 * alpha_val**2) * rng.randn(sample_size)
                  + alpha_val * (X1 + X2))
            X = np.c_[X1, X2, X3]
            y = np.dot(X, coef) + rng.randn(sample_size)

            # get the lasso's coefficients
            _, _, lars_coefs = lars_path(X, y, method='lasso')

            # get the non-negative garotte's coefficients
            ng_coefs, _ = non_negative_garotte_path(X, y, n_alphas=3)

            # test if either model's solution path matches the original model
            if np.any(np.all(ng_coefs.astype(np.bool) == coef.astype(np.bool)[:, np.newaxis], axis=0)):
                ng_path_correct = ng_path_correct + 1.0
            if np.any(np.all(lars_coefs.astype(np.bool) == coef.astype(np.bool)[:, np.newaxis], axis=0)):
                lars_path_correct = lars_path_correct + 1.0

        hits_pers_lars = lars_path_correct/100
        hits_lars.append(hits_pers_lars)
        lars_path_correct = 0
        hits_pers_ng = ng_path_correct/100
        hits_ng.append(hits_pers_ng)
        ng_path_correct = 0

    pl.plot(xrange(25, max_samples, 25), hits_lars, 'r-')
    pl.plot(xrange(25, max_samples, 25), hits_ng, 'b-')
    pl.xlim([0, max_samples])
    pl.ylim([0, 1.1])
    hits_lars = []
    hits_ng = []

pl.show()
