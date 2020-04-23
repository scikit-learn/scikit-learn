import numpy as np

from sklearn.cluster import KMeans
from sklearn.cluster._kmeans import _k_init
from sklearn.utils.extmath import row_norms

from .common import Benchmark, Estimator, Predictor, Transformer
from .datasets import _blobs_dataset, _20newsgroups_highdim_dataset
from .utils import neg_mean_inertia


class KMeans_bench(Predictor, Transformer, Estimator, Benchmark):
    """
    Benchmarks for KMeans.
    """

    param_names = ['representation', 'algorithm', 'n_jobs']
    params = (['dense', 'sparse'], ['full', 'elkan'])

    def setup_cache(self):
        super().setup_cache()

    def setup_cache_(self, params):
        representation, algorithm = params

        if Benchmark.data_size == 'large':
            if representation == 'sparse':
                data = _20newsgroups_highdim_dataset(ngrams=(1, 2))
                n_clusters = 20
            else:
                data = _blobs_dataset()
                n_clusters = 1024
        else:
            if representation == 'sparse':
                data = _20newsgroups_highdim_dataset()
                n_clusters = 20
            else:
                data = _blobs_dataset()
                n_clusters = 256

        estimator = KMeans(n_clusters=n_clusters,
                           algorithm=algorithm,
                           n_init=1,
                           init='random',
                           max_iter=30,
                           tol=1e-16,
                           random_state=0)

        return data, estimator

    def make_scorers(self):
        self.train_scorer = (
            lambda _, __: neg_mean_inertia(self.X,
                                           self.estimator.predict(self.X),
                                           self.estimator.cluster_centers_))
        self.test_scorer = (
            lambda _, __: neg_mean_inertia(self.X_val,
                                           self.estimator.predict(self.X_val),
                                           self.estimator.cluster_centers_))


class KMeansPlusPlus_bench(Benchmark):
    """
    Benchmarks for k-means++ init.
    """

    param_names = ['representation']
    params = (['dense', 'sparse'],)

    def setup(self, *params):
        representation, = params

        if representation == 'sparse':
            data = _20newsgroups_highdim_dataset(ngrams=(1, 2))
            self.n_clusters = 20
        else:
            data = _blobs_dataset()
            self.n_clusters = 256
        self.X, self.X_val, self.y, self.y_val = data

        self.x_squared_norms = row_norms(self.X, squared=True)

    def time_kmeansplusplus(self, *args):
        _k_init(self.X, self.n_clusters, self.x_squared_norms,
                random_state=np.random.RandomState(0))

    def peakmem_kmeansplusplus(self, *args):
        _k_init(self.X, self.n_clusters, self.x_squared_norms,
                random_state=np.random.RandomState(0))
