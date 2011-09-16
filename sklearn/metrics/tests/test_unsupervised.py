import numpy as np
from numpy import linalg
from numpy.testing import assert_true

from ... import datasets
from ..unsupervised import silhouette_coefficient
from ..pairwise import pairwise_distance

def test_silhouette():
    """ Tests the Silhouette Coefficient. """
    dataset = datasets.load_iris()
    X = dataset.data
    y = dataset.target
    D = pairwise_distance(X, metric='euclidean')
    # Given that the actual labels are used, we can assume that S would be
    # positive.
    silhouette = silhouette_coefficient(y, D)
    assert_true(silhouette > 0)
    raise ValueError("Got to here. silhouette")
