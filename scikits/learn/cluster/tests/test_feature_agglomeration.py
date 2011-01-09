"""
Several basic tests for hierarchical clustering procedures

Author : Vincent Michel, 2010
"""

import numpy as np
from scikits.learn.feature_agglomeration import WardAgglomeration
from scikits.learn.feature_extraction.image import img_to_graph


def test_ward_agglomeration():
    """
    Check that we obtain the correct solution in a simplistic case
    """
    np.random.seed(0)
    mask = np.ones([10, 10], dtype=np.bool)
    X = np.random.randn(50, 100)
    adjacency_matrix = img_to_graph(mask, mask)
    feat_agglo = WardAgglomeration(n_clusters=5)
    feat_agglo.fit(X, adjacency_matrix)
    assert(np.size(np.unique(feat_agglo.clustering.label_)) == 5)


if __name__ == '__main__':
    import nose
    nose.run(argv=['', __file__])
