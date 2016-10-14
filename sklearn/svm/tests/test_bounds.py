import warnings

import numpy as np
from scipy import sparse as sp

from sklearn.svm.bounds import l1_min_c
from sklearn.svm import LinearSVC
from sklearn.linear_model.logistic import LogisticRegression

from sklearn.utils.testing import assert_true, raises
from sklearn.utils.testing import assert_raise_message


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

    # loss='l2' should raise ValueError
    assert_raise_message(ValueError, "loss type not in",
                         l1_min_c, dense_X, Y1, "l2")


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


@raises(ValueError)
def test_ill_posed_min_c():
    X = [[0, 0], [0, 0]]
    y = [0, 1]
    l1_min_c(X, y)


@raises(ValueError)
def test_unsupported_loss():
    l1_min_c(dense_X, Y1, 'l1')
