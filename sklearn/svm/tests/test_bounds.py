import nose
from nose.tools import assert_equal, assert_true
from sklearn.utils.testing import clean_warning_registry
import warnings

import numpy as np
from scipy import sparse as sp

from sklearn.svm.bounds import l1_min_c
from sklearn.svm import LinearSVC
from sklearn.linear_model.logistic import LogisticRegression


dense_X = [[-1, 0], [0, 1], [1, 1], [1, 1]]
sparse_X = sp.csr_matrix(dense_X)

Y1 = [0, 1, 1, 1]
Y2 = [2, 1, 0, 0]


def test_l1_min_c():
    losses = ['squared_hinge', 'log']
    Xs = {'sparse': sparse_X, 'dense': dense_X}
    Ys = {'two-classes': Y1, 'multi-class': Y2}
    intercepts = {'no-intercept': {'fit_intercept': False},
                  'fit-intercept': {'fit_intercept': True,
                                    'intercept_scaling': 10}}

    for loss in losses:
        for X_label, X in Xs.items():
            for Y_label, Y in Ys.items():
                for intercept_label, intercept_params in intercepts.items():
                    check = lambda: check_l1_min_c(X, Y, loss,
                                                   **intercept_params)
                    check.description = ('Test l1_min_c loss=%r %s %s %s' %
                                         (loss, X_label, Y_label,
                                          intercept_label))
                    yield check


def test_l2_deprecation():
    clean_warning_registry()
    with warnings.catch_warnings(record=True) as w:
        assert_equal(l1_min_c(dense_X, Y1, "l2"),
                     l1_min_c(dense_X, Y1, "squared_hinge"))
        assert_equal(w[0].category, DeprecationWarning)


def check_l1_min_c(X, y, loss, fit_intercept=True, intercept_scaling=None):
    min_c = l1_min_c(X, y, loss, fit_intercept, intercept_scaling)

    clf = {
        'log': LogisticRegression(penalty='l1'),
        'squared_hinge': LinearSVC(loss='squared_hinge',
                                   penalty='l1', dual=False),
    }[loss]

    clf.fit_intercept = fit_intercept
    clf.intercept_scaling = intercept_scaling

    clf.C = min_c
    clf.fit(X, y)
    assert_true((np.asarray(clf.coef_) == 0).all())
    assert_true((np.asarray(clf.intercept_) == 0).all())

    clf.C = min_c * 1.01
    clf.fit(X, y)
    assert_true((np.asarray(clf.coef_) != 0).any() or
                (np.asarray(clf.intercept_) != 0).any())


@nose.tools.raises(ValueError)
def test_ill_posed_min_c():
    X = [[0, 0], [0, 0]]
    y = [0, 1]
    l1_min_c(X, y)


@nose.tools.raises(ValueError)
def test_unsupported_loss():
    l1_min_c(dense_X, Y1, 'l1')
