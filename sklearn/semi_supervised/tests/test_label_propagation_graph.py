""" test the label propagation module """

import nose
import numpy as np

from sklearn.semi_supervised import label_propagation
from scipy.sparse import dok_matrix
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal


ESTIMATORS = [
    (label_propagation.LabelPropagation, {'kernel': None}),
    (label_propagation.LabelSpreading, {'kernel': None}),
]

graph = [[0, 1, 0], [1, 0, 0], [0, 0, 1]]
labels = [0, -1, 1]

def test_fit_transduction():
    for estimator, parameters in ESTIMATORS:
        clf = estimator(**parameters).fit(graph, labels)
        nose.tools.assert_equal(clf.transduction_[2], 1)

def test_fit_transduction_sparse():
    graph_sparse = dok_matrix(graph)
    for estimator, parameters in ESTIMATORS:
        clf = estimator(**parameters).fit(graph_sparse, labels)
        nose.tools.assert_equal(clf.transduction_[2], 1)

def test_distribution():
    for estimator, parameters in ESTIMATORS:
        clf = estimator(**parameters).fit(graph, labels)
        assert_array_almost_equal(np.asarray(clf.label_distributions_[1]),
                                      np.array([1, 0]), 2)


def test_predict():
    for estimator, parameters in ESTIMATORS:
        clf = estimator(**parameters).fit(graph, labels)
        assert_array_equal(clf.predict([[1, 0, 0]]), np.array([0, 0, 0]))


def test_predict_proba():
    for estimator, parameters in ESTIMATORS:
        clf = estimator(**parameters).fit(graph, labels)
        assert_array_almost_equal(clf.predict_proba([[1, 0, 0]]),
                                  np.array([[0.5, 0.5]]))