""" test the label propagation module """
from .. import label_propagation
from numpy.testing import assert_array_almost_equal
from nose.tools import assert_equal, assert_raises
from numpy.testing import assert_array_equal

samples = [[1,0],[0,1],[1,1]]
labels = [[0,1],[1,0], None]

def test_base_multiclass_fit():
    lp = label_propagation.LabelPropagation()
    lp.fit(samples, labels)
    assert lp.y[2] == [.6, .5]

def test_base_binary_fit():
    samples = [[1,0],[0,1],[1,1]]
    labels = [-1, 1, None]

    lp = label_propagation.LabelPropagation()
    lp.fit(samples, labels)
    assert lp.y[2] == 0

#def test_predict_probabilities(self):
#
#def test_predict(self)
#
#def test_label_propagation(self):
#def test_label_spreading(self):

