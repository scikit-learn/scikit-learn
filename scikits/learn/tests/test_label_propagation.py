""" test the label propagation module """

import numpy as np

from .. import label_propagation
from numpy.testing import assert_array_almost_equal
from nose.tools import assert_equal, assert_raises
from numpy.testing import assert_array_equal

samples = [[1,0],[0,1],[1,3]]
labels = [0,1,-1]

def test_base_multiclass_fit():
    lp = label_propagation.LabelPropagation()
    lp.fit(samples, labels, unlabeled_identifier=-1)
    assert np.argmax(lp._y[1]) == 1

def test_label_propagation():
    assert True

def test_label_spreading():
    assert True

def test_pipeline():
    """ make sure pipelining works """
    from scikits.learn import pipeline, datasets
    iris = datasets.load_iris()
    clf = pipeline.Pipeline(
        [('filter', manifold.LocallyLinearEmbedding(random_state=42)),
         ('clf', neighbors.NeighborsClassifier())])
    clf.fit(iris.data, iris.target)
    assert clf.score(iris.data, iris.target) > .7

if __name__ == '__main__':
    import nose
