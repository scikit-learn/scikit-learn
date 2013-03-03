import itertools

import numpy as np
from numpy.testing import assert_allclose

from scipy.spatial.distance import cdist, pdist, squareform
from sklearn.neighbors.dist_metrics import DistanceMetric


class TestMetrics:
    def __init__(self, n1=20, n2=25, d=4, zero_frac=0.5,
                 rseed=0, dtype=np.float64):
        np.random.seed(rseed)
        self.X1 = np.random.random((n1, d)).astype(dtype)
        self.X2 = np.random.random((n2, d)).astype(dtype)

        # make boolean arrays: ones and zeros
        self.X1_bool = self.X1.round(0)
        self.X2_bool = self.X2.round(0)

        V = np.random.random((d, d))
        VI = np.dot(V, V.T)

        self.metrics = {'euclidean':{},
                        'cityblock':{},
                        'minkowski':dict(p=(1, 1.5, 2, 3)),
                        'chebyshev':{},
                        'seuclidean':dict(V=(np.random.random(d),)),
                        'wminkowski':dict(p=(1, 1.5, 3),
                                          w=(np.random.random(d),)),
                        'mahalanobis':dict(VI=(VI,)),
                        'hamming':{},
                        'canberra':{},
                        'braycurtis':{}}

        self.bool_metrics = ['matching', 'jaccard', 'dice',
                             'kulsinski', 'rogerstanimoto', 'russellrao',
                             'sokalmichener', 'sokalsneath']

    def test_cdist(self):
        for metric, argdict in self.metrics.iteritems():
            keys = argdict.keys()
            for vals in itertools.product(*argdict.values()):
                kwargs = dict(zip(keys, vals))
                D_true = cdist(self.X1, self.X2, metric, **kwargs)
                yield self.check_cdist, metric, kwargs, D_true

        for metric in self.bool_metrics:
            D_true = cdist(self.X1_bool, self.X2_bool, metric)
            yield self.check_cdist_bool, metric, D_true
            
    def check_cdist(self, metric, kwargs, D_true):
        dm = DistanceMetric.get_metric(metric, **kwargs)
        D12 = dm.pairwise(self.X1, self.X2)
        assert_allclose(D12, D_true)

    def check_cdist_bool(self, metric, D_true):
        dm = DistanceMetric.get_metric(metric)
        D12 = dm.pairwise(self.X1_bool, self.X2_bool)
        assert_allclose(D12, D_true)

    def test_pdist(self):
        for metric, argdict in self.metrics.iteritems():
            keys = argdict.keys()
            for vals in itertools.product(*argdict.values()):
                kwargs = dict(zip(keys, vals))
                D_true = pdist(self.X1, metric, **kwargs)
                Dsq_true = squareform(D_true)
                yield self.check_pdist, metric, kwargs, Dsq_true

    def check_pdist(self, metric, kwargs, D_true):
        dm = DistanceMetric.get_metric(metric, **kwargs)
        D12 = dm.pairwise(self.X1)
        assert_allclose(D12, D_true)        

        
if __name__ == '__main__':
    import nose
    nose.runmodule()
