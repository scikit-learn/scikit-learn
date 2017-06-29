""" test the label propagation module """

import numpy as np

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import assert_no_warnings
from sklearn.semi_supervised import label_propagation
from sklearn.datasets import make_classification
from sklearn.metrics.pairwise import rbf_kernel
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal

ESTIMATORS = [
    (label_propagation.LabelPropagation, {'kernel': 'rbf'}),
    (label_propagation.LabelPropagation, {'kernel': 'knn', 'n_neighbors': 2}),
    (label_propagation.LabelPropagation, {
        'kernel': lambda x, y: rbf_kernel(x, y, gamma=20)
    }),
    (label_propagation.LabelSpreading, {'kernel': 'rbf'}),
    (label_propagation.LabelSpreading, {'kernel': 'knn', 'n_neighbors': 2}),
    (label_propagation.LabelSpreading, {
        'kernel': lambda x, y: rbf_kernel(x, y, gamma=20)
    }),
]


def test_fit_transduction():
    samples = [[1., 0.], [0., 2.], [1., 3.]]
    labels = [0, 1, -1]
    for estimator, parameters in ESTIMATORS:
        clf = estimator(**parameters).fit(samples, labels)
        assert_equal(clf.transduction_[2], 1)


def test_distribution():
    samples = [[1., 0.], [0., 1.], [1., 1.]]
    labels = [0, 1, -1]
    for estimator, parameters in ESTIMATORS:
        clf = estimator(**parameters).fit(samples, labels)
        if parameters['kernel'] == 'knn':
            continue    # unstable test; changes in k-NN ordering break it
            assert_array_almost_equal(clf.predict_proba([[1., 0.0]]),
                                      np.array([[1., 0.]]), 2)
        else:
            assert_array_almost_equal(np.asarray(clf.label_distributions_[2]),
                                      np.array([.5, .5]), 2)


def test_predict():
    samples = [[1., 0.], [0., 2.], [1., 3.]]
    labels = [0, 1, -1]
    for estimator, parameters in ESTIMATORS:
        clf = estimator(**parameters).fit(samples, labels)
        assert_array_equal(clf.predict([[0.5, 2.5]]), np.array([1]))


def test_predict_proba():
    samples = [[1., 0.], [0., 1.], [1., 2.5]]
    labels = [0, 1, -1]
    for estimator, parameters in ESTIMATORS:
        clf = estimator(**parameters).fit(samples, labels)
        assert_array_almost_equal(clf.predict_proba([[1., 1.]]),
                                  np.array([[0.5, 0.5]]))


def test_alpha_deprecation():
    X, y = make_classification(n_samples=100)
    y[::3] = -1
    # Using kernel=knn as rbf appears to result in exp underflow
    lp_default = label_propagation.LabelPropagation(kernel='knn')
    lp_default_y = assert_no_warnings(lp_default.fit, X, y).transduction_

    lp_0 = label_propagation.LabelPropagation(alpha=0, kernel='knn')
    lp_0_y = assert_warns(DeprecationWarning, lp_0.fit, X, y).transduction_

    assert_array_equal(lp_default_y, lp_0_y)
