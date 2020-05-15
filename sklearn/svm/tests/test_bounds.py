import numpy as np
from scipy import sparse as sp
from scipy import stats

import pytest

from sklearn.svm._bounds import l1_min_c
from sklearn.svm import LinearSVC
from sklearn.linear_model import LogisticRegression
from sklearn.svm._newrand import bounded_rand_int_wrap

from sklearn.utils._testing import assert_raise_message


dense_X = [[-1, 0], [0, 1], [1, 1], [1, 1]]
sparse_X = sp.csr_matrix(dense_X)

Y1 = [0, 1, 1, 1]
Y2 = [2, 1, 0, 0]


@pytest.mark.parametrize('loss', ['squared_hinge', 'log'])
@pytest.mark.parametrize('X_label', ['sparse', 'dense'])
@pytest.mark.parametrize('Y_label', ['two-classes', 'multi-class'])
@pytest.mark.parametrize('intercept_label', ['no-intercept', 'fit-intercept'])
def test_l1_min_c(loss, X_label, Y_label, intercept_label):
    Xs = {'sparse': sparse_X, 'dense': dense_X}
    Ys = {'two-classes': Y1, 'multi-class': Y2}
    intercepts = {'no-intercept': {'fit_intercept': False},
                  'fit-intercept': {'fit_intercept': True,
                                    'intercept_scaling': 10}}

    X = Xs[X_label]
    Y = Ys[Y_label]
    intercept_params = intercepts[intercept_label]
    check_l1_min_c(X, Y, loss, **intercept_params)


def test_l1_min_c_l2_loss():
    # loss='l2' should raise ValueError
    assert_raise_message(ValueError, "loss type not in",
                         l1_min_c, dense_X, Y1, loss="l2")


def check_l1_min_c(X, y, loss, fit_intercept=True, intercept_scaling=None):
    min_c = l1_min_c(X, y, loss=loss, fit_intercept=fit_intercept,
                     intercept_scaling=intercept_scaling)

    clf = {
        'log': LogisticRegression(penalty='l1', solver='liblinear'),
        'squared_hinge': LinearSVC(loss='squared_hinge',
                                   penalty='l1', dual=False),
    }[loss]

    clf.fit_intercept = fit_intercept
    clf.intercept_scaling = intercept_scaling

    clf.C = min_c
    clf.fit(X, y)
    assert (np.asarray(clf.coef_) == 0).all()
    assert (np.asarray(clf.intercept_) == 0).all()

    clf.C = min_c * 1.01
    clf.fit(X, y)
    assert ((np.asarray(clf.coef_) != 0).any() or
            (np.asarray(clf.intercept_) != 0).any())


def test_ill_posed_min_c():
    X = [[0, 0], [0, 0]]
    y = [0, 1]
    with pytest.raises(ValueError):
        l1_min_c(X, y)


def test_unsupported_loss():
    with pytest.raises(ValueError):
        l1_min_c(dense_X, Y1, loss='l1')


# TODO test deterministic results on differing systems
# def test_set_seed():
#     def _test(seed, val):
#         # TODO currently causes 99% coverage in test_bounds.py
#         if seed is not None:
#             set_seed_wrap(seed)
#         x = bounded_rand_int_wrap(100)
#         assert(x == val), 'Expected {} but got {} instead'.format(val, x)

#     # TODO should be default seed for std::mt19937,
#     # but different results on jupyter lab
#     # _test(None, 81)
#     _test(5489, 81)  # default seed for std::mt19937
#     _test(0, 54)
#     _test(4294967295, 9)  # max unsigned int size


_MAX_INT = 2147483647
@pytest.mark.parametrize('orig_range, n_pts',
                         [(_MAX_INT, 10000), (100, 10)])
def test_bounded_rand_int(orig_range, n_pts):
    n_iter = 100
    ks_pvals = []
    uniform_dist = stats.uniform(loc=0, scale=orig_range)
    # perform multiple samplings to make chance of outlier sampling negligible
    for _ in range(n_iter):
        # Deterministic random sampling
        sample = [bounded_rand_int_wrap(orig_range) for _ in range(n_pts)]
        res = stats.kstest(sample, uniform_dist.cdf)
        ks_pvals.append(res.pvalue)
    min_10pct_pval = np.quantile(ks_pvals, q=0.1)
    # lower 10th quantile pvalue <= 0.05 means that the test rejects the
    # null hypothesis that the sample came from the uniform distribution
    assert(min_10pct_pval > 0.05),\
        'lower 10th quantile p-value of {} not > 0.05, '\
        'therefore distribution isn\'t sufficiently uniform'\
        .format(min_10pct_pval)
