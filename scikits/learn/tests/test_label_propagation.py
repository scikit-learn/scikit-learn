""" test the label propagation module """

import numpy as np

from .. import label_propagation
from numpy.testing import assert_array_almost_equal
from nose.tools import assert_equal, assert_raises
from numpy.testing import assert_array_equal


def test_label_propagation_fit():
    samples = [[1,0],[0,1],[1,3]]
    labels = [0,1,-1]
    lp = label_propagation.LabelPropagation()
    lp.fit(samples, labels, unlabeled_identifier=-1)
    assert lp.transduction[2] == 1

def test_label_spreading_fit():
    samples = [[1,0],[0,1],[1,3]]
    labels = [0,1,-1]
    lp = label_propagation.LabelSpreading()
    lp.fit(samples, labels, unlabeled_identifier=-1)
    assert lp.transduction[2] == 1

def test_label_prop():
    X_, y_ = test_dataset_classif(n_samples=200, n_features=100, seed=0)
    unlabeled = np.random.randint(0, 200

def test_string_labels():
    samples = [[1,0],[0,1],[1,3]]
    labels = ['banana', 'orange', 'unabeled']
    lp = label_propagation.LabelPropagation()
    lp.fit(samples, labels, unlabeled_identifier='unlabeled')
    assert lp.transduction[2] == 'unlabeled'

def test_distribution():
    samples = [[1,0],[0,1],[1,3]]
    labels = [0,1,-1]
    lp = label_propagation.LabelPropagation()
    lp.fit(samples, labels, unlabeled_identifier=-1)
    assert_array_equal(lp._y[2], [.5,.5]

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
