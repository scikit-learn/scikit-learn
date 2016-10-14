""" test the label propagation module """

import numpy as np

from sklearn.utils.testing import assert_equal
from sklearn.semi_supervised import label_propagation
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal


ESTIMATORS = [
    (label_propagation.LabelPropagation, {'kernel': 'rbf'}),
    (label_propagation.LabelPropagation, {'kernel': 'knn', 'n_neighbors': 2}),
    (label_propagation.LabelSpreading, {'kernel': 'rbf'}),
    (label_propagation.LabelSpreading, {'kernel': 'knn', 'n_neighbors': 2})
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
