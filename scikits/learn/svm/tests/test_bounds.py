import numpy as np
from scipy import sparse as sp

from numpy.testing import assert_array_equal, \
     assert_array_almost_equal, assert_almost_equal

from scikits.learn.svm.bounds import l1_min_c
from scikits.learn.svm import LinearSVC
from scikits.learn.linear_model.logistic import LogisticRegression
from scikits.learn.svm.sparse import LinearSVC as SparseLinearSVC
from scikits.learn.linear_model.sparse.logistic import LogisticRegression as \
                                                       SparseLogisticRegression


dense_X = [[-1, 0], [0, 1], [1, 1], [1, 1]]
sparse_X = sp.csr_matrix(dense_X)

Y1 = [0, 1, 1, 1]
Y2 = [2, 1, 0, 0]


def test_l1_min_c():
    losses = ['l2', 'log']
    Xs = {'sparse': sparse_X, 'dense': dense_X}
    Ys = {'two-classes': Y1, 'multi-class': Y2}
    intercepts = {'no-intercept':  {'fit_intercept': False},
                  'fit-intercept': {'fit_intercept': True,
                                    'intercept_scaling': 10}}

    for loss in losses:
        for X_label, X in Xs.items():
            for Y_label, Y in Ys.items():
                for intercept_label, intercept_params in intercepts.items():
                    desc = 'l1_min_c loss=%r %s %s %s' % \
                                      (loss, X_label, Y_label, intercept_label)
                    check = lambda X, Y, loss, intercept_params: \
                                   check_l1_min_c(X, Y, loss, intercept_params)
                    check.description = desc
                    yield check, X, Y, loss, intercept_params


def check_l1_min_c(X, y, loss, intercept_params):
    min_c = l1_min_c(X, y, loss, **intercept_params)

    if loss == 'l2':
        if sp.issparse(X):
            algo = SparseLinearSVC
        else:
            algo = LinearSVC
        clf = algo(loss='l2', penalty='l1', dual=False, **intercept_params)
    else:
        if sp.issparse(X):
            algo = SparseLogisticRegression
        else:
            algo = LogisticRegression
        clf = algo(penalty='l1', **intercept_params)

    clf.C = min_c
    clf.fit(X, y)
    assert (np.asanyarray(clf.coef_) == 0).all()
    assert (np.asanyarray(clf.intercept_) == 0).all()

    clf.C = min_c * 1.01
    clf.fit(X, y)
    assert (np.asanyarray(clf.coef_) != 0).any() or \
           (np.asanyarray(clf.intercept_) != 0).any()
