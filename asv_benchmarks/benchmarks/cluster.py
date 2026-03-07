from sklearn.cluster import HDBSCAN, KMeans, MiniBatchKMeans
from sklearn.datasets import make_blobs

from .common import Benchmark, Estimator, Predictor, Transformer
from .datasets import _20newsgroups_highdim_dataset, _blobs_dataset
from .utils import neg_mean_inertia


class KMeansBenchmark(Predictor, Transformer, Estimator, Benchmark):
    """
    Benchmarks for KMeans.
    """

    param_names = ["representation", "algorithm", "init"]
    params = (["dense", "sparse"], ["lloyd", "elkan"], ["random", "k-means++"])

    def setup_cache(self):
        super().setup_cache()

    def make_data(self, params):
        representation, algorithm, init = params

        if representation == "sparse":
            data = _20newsgroups_highdim_dataset(n_samples=8000)
        else:
            data = _blobs_dataset(n_clusters=20)

        return data

    def make_estimator(self, params):
        representation, algorithm, init = params

        max_iter = 30 if representation == "sparse" else 100

        estimator = KMeans(
            n_clusters=20,
            algorithm=algorithm,
            init=init,
            n_init=1,
            max_iter=max_iter,
            tol=0,
            random_state=0,
        )

        return estimator

    def make_scorers(self):
        self.train_scorer = lambda _, __: neg_mean_inertia(
            self.X, self.estimator.predict(self.X), self.estimator.cluster_centers_
        )
        self.test_scorer = lambda _, __: neg_mean_inertia(
            self.X_val,
            self.estimator.predict(self.X_val),
            self.estimator.cluster_centers_,
        )


class MiniBatchKMeansBenchmark(Predictor, Transformer, Estimator, Benchmark):
    """
    Benchmarks for MiniBatchKMeans.
    """

    param_names = ["representation", "init"]
    params = (["dense", "sparse"], ["random", "k-means++"])

    def setup_cache(self):
        super().setup_cache()

    def make_data(self, params):
        representation, init = params

        if representation == "sparse":
            data = _20newsgroups_highdim_dataset()
        else:
            data = _blobs_dataset(n_clusters=20)

        return data

    def make_estimator(self, params):
        representation, init = params

        max_iter = 5 if representation == "sparse" else 2

        estimator = MiniBatchKMeans(
            n_clusters=20,
            init=init,
            n_init=1,
            max_iter=max_iter,
            batch_size=1000,
            max_no_improvement=None,
            compute_labels=False,
            random_state=0,
        )

        return estimator

    def make_scorers(self):
        self.train_scorer = lambda _, __: neg_mean_inertia(
            self.X, self.estimator.predict(self.X), self.estimator.cluster_centers_
        )
        self.test_scorer = lambda _, __: neg_mean_inertia(
            self.X_val,
            self.estimator.predict(self.X_val),
            self.estimator.cluster_centers_,
        )


class HDBSCANMSTBenchmark:
    """Benchmark Prim's and Boruvka MST backends in HDBSCAN."""

    param_names = ["mst_algorithm"]
    params = (["prims", "boruvka_exact"],)
    timeout = 180.0

    def setup(self, mst_algorithm):
        self.X, _ = make_blobs(
            n_samples=10000,
            n_features=20,
            centers=5,
            cluster_std=0.6,
            random_state=0,
        )
        self.estimator = HDBSCAN(
            algorithm="auto",
            mst_algorithm=mst_algorithm,
            n_jobs=-1,
            copy=False,
        )

    def time_fit(self, mst_algorithm):
        self.estimator.fit(self.X)
