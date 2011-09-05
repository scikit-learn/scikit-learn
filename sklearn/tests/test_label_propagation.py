""" test the label propagation module """

import numpy as np

from .. import label_propagation
from numpy.testing import assert_array_almost_equal
from nose.tools import assert_equal, assert_raises
from numpy.testing import assert_array_equal
from numpy.testing import assert_equal

from StringIO import StringIO


def test_label_propagation_fit():
    samples = [[1., 0.], [0., 1.], [1., 3.]]
    labels = [0, 1, -1]
    lp = label_propagation.LabelPropagation()
    lp.fit(samples, labels)
    assert lp.transduction_[2] == 1


def test_label_spreading_fit():
    samples = [[1., 0.], [0., 1.], [1., 3.]]
    labels = [0, 1, -1]
    lp = label_propagation.LabelSpreading()
    lp.fit(samples, labels)
    assert lp.transduction_[2] == 1


def test_string_labels():
    samples = [[1., 0.], [0., 1.], [1., 3.]]
    labels = ['banana', 'orange', 'unlabeled']
    lp = label_propagation.LabelPropagation(unlabeled_identifier='unlabeled')
    lp.fit(samples, labels)
    assert lp.transduction_[2] == 'orange'


def test_distribution():
    samples = [[1., 0.], [0., 1.], [1., 1.]]
    labels = [0, 1, -1]
    lp = label_propagation.LabelPropagation()
    lp.fit(samples, labels)
    assert_array_almost_equal(np.asarray(lp.y_[2]),
             np.array([  2.06115361e-09,   4.12230722e-09]))


def test_predict():
    samples = [[1., 0.], [0., 1.], [1., 3.]]
    labels = [0, 1, -1]
    lp = label_propagation.LabelPropagation()
    lp.fit(samples, labels)
    assert lp.predict([[1,1]]) == np.array([1])


def test_predict_proba():
    samples = [[1., 0.], [0., 1.], [1., 3.]]
    labels = [0, 1, -1]
    lp = label_propagation.LabelPropagation()
    lp.fit(samples, labels)
    assert_array_almost_equal(lp.predict_proba([[1,1]]),
            np.array([[  8.75651076e-27,   2.06115362e-09]]))


if __name__ == '__main__':
    import nose
