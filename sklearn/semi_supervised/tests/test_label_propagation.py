""" test the label propagation module """

import numpy as np
import pytest

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_no_warnings
from sklearn.semi_supervised import label_propagation
from sklearn.metrics.pairwise import rbf_kernel
from sklearn.datasets import make_classification
from sklearn.exceptions import ConvergenceWarning
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal


ESTIMATORS = (label_propagation.LabelPropagation,
              label_propagation.LabelSpreading)
ESTIMATOR_PARAMS = ({'kernel': 'rbf'},
                    {'kernel': 'knn', 'n_neighbors': 2},
                    {'kernel': lambda x, y: rbf_kernel(x, y, gamma=20)})


@pytest.mark.parametrize('estimator', ESTIMATORS)
@pytest.mark.parametrize('parameters', ESTIMATOR_PARAMS)
def test_fit_transduction(estimator, parameters):
    samples = [[1., 0.], [0., 2.], [1., 3.]]
    labels = [0, 1, -1]
    clf = estimator(**parameters).fit(samples, labels)
    assert_equal(clf.transduction_[2], 1)


@pytest.mark.parametrize('estimator', ESTIMATORS)
@pytest.mark.parametrize('parameters', ESTIMATOR_PARAMS)
def test_distribution(estimator, parameters):
    samples = [[1., 0.], [0., 1.], [1., 1.]]
    labels = [0, 1, -1]
    clf = estimator(**parameters).fit(samples, labels)
    if parameters['kernel'] == 'knn':
        pytest.skip('Unstable test; changes in k-NN ordering break it')
        assert_array_almost_equal(clf.predict_proba([[1., 0.0]]),
                                  np.array([[1., 0.]]), 2)
    else:
        assert_array_almost_equal(np.asarray(clf.label_distributions_[2]),
                                  np.array([.5, .5]), 2)


@pytest.mark.parametrize('estimator', ESTIMATORS)
@pytest.mark.parametrize('parameters', ESTIMATOR_PARAMS)
def test_predict(estimator, parameters):
    samples = [[1., 0.], [0., 2.], [1., 3.]]
    labels = [0, 1, -1]
    clf = estimator(**parameters).fit(samples, labels)
    assert_array_equal(clf.predict([[0.5, 2.5]]), np.array([1]))


@pytest.mark.parametrize('estimator', ESTIMATORS)
@pytest.mark.parametrize('parameters', ESTIMATOR_PARAMS)
def test_predict_proba(estimator, parameters):
    samples = [[1., 0.], [0., 1.], [1., 2.5]]
    labels = [0, 1, -1]
    clf = estimator(**parameters).fit(samples, labels)
    assert_array_almost_equal(clf.predict_proba([[1., 1.]]),
                              np.array([[0.5, 0.5]]))


def test_alpha_deprecation():
    X, y = make_classification(n_samples=100)
    y[::3] = -1

    lp_default = label_propagation.LabelPropagation(kernel='rbf', gamma=0.1)
    lp_default_y = lp_default.fit(X, y).transduction_

    lp_0 = label_propagation.LabelPropagation(alpha=0, kernel='rbf', gamma=0.1)
    lp_0_y = assert_warns(DeprecationWarning, lp_0.fit, X, y).transduction_

    assert_array_equal(lp_default_y, lp_0_y)


def test_label_spreading_closed_form():
    n_classes = 2
    X, y = make_classification(n_classes=n_classes, n_samples=200,
                               random_state=0)
    y[::3] = -1
    clf = label_propagation.LabelSpreading().fit(X, y)
    # adopting notation from Zhou et al (2004):
    S = clf._build_graph()
    Y = np.zeros((len(y), n_classes + 1))
    Y[np.arange(len(y)), y] = 1
    Y = Y[:, :-1]
    for alpha in [0.1, 0.3, 0.5, 0.7, 0.9]:
        expected = np.dot(np.linalg.inv(np.eye(len(S)) - alpha * S), Y)
        expected /= expected.sum(axis=1)[:, np.newaxis]
        clf = label_propagation.LabelSpreading(max_iter=10000, alpha=alpha)
        clf.fit(X, y)
        assert_array_almost_equal(expected, clf.label_distributions_, 4)


def test_label_propagation_closed_form():
    n_classes = 2
    X, y = make_classification(n_classes=n_classes, n_samples=200,
                               random_state=0)
    y[::3] = -1
    Y = np.zeros((len(y), n_classes + 1))
    Y[np.arange(len(y)), y] = 1
    unlabelled_idx = Y[:, (-1,)].nonzero()[0]
    labelled_idx = (Y[:, (-1,)] == 0).nonzero()[0]

    clf = label_propagation.LabelPropagation(max_iter=10000,
                                             gamma=0.1)
    clf.fit(X, y)
    # adopting notation from Zhu et al 2002
    T_bar = clf._build_graph()
    Tuu = T_bar[np.meshgrid(unlabelled_idx, unlabelled_idx, indexing='ij')]
    Tul = T_bar[np.meshgrid(unlabelled_idx, labelled_idx, indexing='ij')]
    Y = Y[:, :-1]
    Y_l = Y[labelled_idx, :]
    Y_u = np.dot(np.dot(np.linalg.inv(np.eye(Tuu.shape[0]) - Tuu), Tul), Y_l)

    expected = Y.copy()
    expected[unlabelled_idx, :] = Y_u
    expected /= expected.sum(axis=1)[:, np.newaxis]

    assert_array_almost_equal(expected, clf.label_distributions_, 4)


def test_valid_alpha():
    n_classes = 2
    X, y = make_classification(n_classes=n_classes, n_samples=200,
                               random_state=0)
    for alpha in [-0.1, 0, 1, 1.1, None]:
        assert_raises(ValueError,
                      lambda **kwargs:
                      label_propagation.LabelSpreading(**kwargs).fit(X, y),
                      alpha=alpha)


def test_convergence_speed():
    # This is a non-regression test for #5774
    X = np.array([[1., 0.], [0., 1.], [1., 2.5]])
    y = np.array([0, 1, -1])
    mdl = label_propagation.LabelSpreading(kernel='rbf', max_iter=5000)
    mdl.fit(X, y)

    # this should converge quickly:
    assert mdl.n_iter_ < 10
    assert_array_equal(mdl.predict(X), [0, 1, 1])


def test_convergence_warning():
    # This is a non-regression test for #5774
    X = np.array([[1., 0.], [0., 1.], [1., 2.5]])
    y = np.array([0, 1, -1])
    mdl = label_propagation.LabelSpreading(kernel='rbf', max_iter=1)
    assert_warns(ConvergenceWarning, mdl.fit, X, y)
    assert_equal(mdl.n_iter_, mdl.max_iter)

    mdl = label_propagation.LabelPropagation(kernel='rbf', max_iter=1)
    assert_warns(ConvergenceWarning, mdl.fit, X, y)
    assert_equal(mdl.n_iter_, mdl.max_iter)

    mdl = label_propagation.LabelSpreading(kernel='rbf', max_iter=500)
    assert_no_warnings(mdl.fit, X, y)

    mdl = label_propagation.LabelPropagation(kernel='rbf', max_iter=500)
    assert_no_warnings(mdl.fit, X, y)
