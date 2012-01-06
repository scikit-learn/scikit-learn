""" test the label propagation module """

from StringIO import StringIO

import numpy as np

from ..label_propagation import LabelPropagation, LabelSpreading
from numpy.testing import assert_array_almost_equal
from nose.tools import assert_equal, assert_raises
from numpy.testing import assert_array_equal
from numpy.testing import assert_equal


def test_label_propagation_fit():
    samples = [[1., 0.], [0., 1.], [1., 3.]]
    labels = [0, 1, -1]
    lp = LabelPropagation(unlabeled_identifier=-1)
    lp.fit(samples, labels)
    assert lp.transduction_[2] == 1


def test_label_spreading_fit():
    samples = [[1., 0.], [0., 1.], [1., 3.]]
    labels = [0, 1, -1]
    lp = LabelSpreading(unlabeled_identifier=-1)
    lp.fit(samples, labels)
    assert lp.transduction_[2] == 1


def test_string_labels():
    samples = [[1., 0.], [0., 1.], [1., 3.]]
    labels = ['banana', 'orange', 'unlabeled']
    lp = LabelPropagation(unlabeled_identifier='unlabeled')
    lp.fit(samples, labels)
    assert lp.transduction_[2] == 'orange'


def test_distribution():
    samples = [[1., 0.], [0., 1.], [1., 1.]]
    labels = [0, 1, -1]
    lp = LabelPropagation(unlabeled_identifier=-1)
    lp.fit(samples, labels)
    assert_array_almost_equal(np.asarray(lp.y_[2]), \
            np.array([[2.06115361e-09, 2.06115361e-09]]))


if __name__ == '__main__':
    import nose
