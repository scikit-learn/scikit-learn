""" test the label propagation module """

import numpy as np

from .. import label_propagation
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal


ESTIMATORS = [
    label_propagation.LabelPropagation,
    label_propagation.LabelSpreading,
]


def test_fit_transduction():
    samples = [[1., 0.], [0., 2.], [1., 3.]]
    labels = [0, 1, -1]
    for estimator in ESTIMATORS:
        clf = estimator().fit(samples, labels)
        assert clf.transduction_[2] == 1


def test_string_labels():
    samples = [[1., 0.], [0., 2.], [1., 3.]]
    labels = ['banana', 'orange', 'unlabeled']
    for estimator in ESTIMATORS:
        clf = estimator(unlabeled_identifier='unlabeled').fit(samples, labels)
        assert clf.transduction_[2] == 'orange'


def test_distribution():
    samples = [[1., 0.], [0., 1.], [1., 1.]]
    labels = [0, 1, -1]
    clf = label_propagation.LabelPropagation().fit(samples, labels)
    print clf.y_
    assert_array_almost_equal(np.asarray(clf.y_[2]), np.zeros(2), 2)

def test_predict():
    samples = [[1., 0.], [0., 2.], [1., 3.]]
    labels = [0, 1, -1]
    for estimator in ESTIMATORS:
        clf = estimator().fit(samples, labels)
        assert_array_equal(clf.predict([[0.5, 2.5]]), np.array([1]))


def test_predict_proba():
    samples = [[1., 0.], [0., 1.], [1., 2.5]]
    labels = [0, 1, -1]

    for estimator in ESTIMATORS:
        clf = estimator().fit(samples, labels)
        assert_array_almost_equal(clf.predict_proba([[1., 1.]]),
                np.array([[0.5, 0.5]]))


if __name__ == '__main__':
    import nose
