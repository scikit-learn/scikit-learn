import numpy as np
from scipy.spatial.distance import cdist
from scipy.sparse import rand as sparse_rand

from .common import Benchmark

from sklearn.metrics._pairwise_distances_reduction import ArgKmin, RadiusNeighbors

# To run benchmarks defined this file, between for instance your <current_branch>
# and upstream/main, use:
#
#     asv continuous -b PairwiseDistancesReductions upstream/main <current_branch>


class PairwiseDistancesReductionsBenchmark(Benchmark):

    param_names = [
        "n_train",
        "n_test",
        "n_features",
        "metric",
        "strategy",
        "dtype",
        "X_train",
        "X_test",
    ]
    params = [
        [1000, 10_000, int(1e7)],
        [1000, 10_000, 100_000],
        [100],
        ["euclidean", "manhattan"],
        ["auto", "parallel_on_X", "parallel_on_Y"],
        [np.float32, np.float64],
        ["dense", "csr"],
        ["dense", "csr"],
    ]

    def setup(
        self, n_train, n_test, n_features, metric, strategy, dtype, X_train, X_test
    ):
        rng = np.random.RandomState(0)
        self.X_train = (
            rng.rand(n_train, n_features).astype(dtype)
            if X_train == "dense"
            else sparse_rand(
                n_train,
                n_features,
                density=0.05,
                format="csr",
                dtype=dtype,
                random_state=rng,
            )
        )
        self.X_test = (
            rng.rand(n_test, n_features).astype(dtype)
            if X_test == "dense"
            else sparse_rand(
                n_test,
                n_features,
                density=0.05,
                format="csr",
                dtype=dtype,
                random_state=rng,
            )
        )

        self.y_train = rng.randint(low=-1, high=1, size=(n_train,))
        self.metric = metric
        self.strategy = strategy

        self.k = 10

        # Motive: radius has to be scaled with the number of feature
        # Defining it as the 0.001-quantile allows to have in expectation
        # a constant amount of neighbors, regardless of the value of n_features.
        dist_mat = cdist(
            (self.X_train if X_train == "dense" else self.X_train.toarray())[:1000],
            (self.X_test if X_test == "dense" else self.X_test.toarray())[:10],
        )

        self.radius = np.quantile(a=dist_mat.ravel(), q=0.001)

    def time_ArgKmin(
        self,
        n_train,
        n_test,
        n_features,
        metric,
        strategy,
        dtype,
        X_train,
        X_test,
    ):

        ArgKmin.compute(
            X=self.X_test,
            Y=self.X_train,
            k=10,
            metric=self.metric,
            return_distance=True,
            strategy=self.strategy,
        )

    def peakmem_ArgKmin(
        self,
        n_train,
        n_test,
        n_features,
        metric,
        strategy,
        dtype,
        X_train,
        X_test,
    ):
        self.time_ArgKmin(
            n_train,
            n_test,
            n_features,
            metric,
            strategy,
            dtype,
            X_train,
            X_test,
        )

    def time_RadiusNeighbors(
        self, n_train, n_test, n_features, metric, strategy, dtype, X_train, X_test
    ):

        RadiusNeighbors.compute(
            X=self.X_test,
            Y=self.X_train,
            radius=self.radius,
            metric=self.metric,
            return_distance=True,
            strategy=self.strategy,
        )

    def peakmem_RadiusNeighbors(
        self, n_train, n_test, n_features, metric, strategy, dtype, X_train, X_test
    ):
        self.time_RadiusNeighbors(
            n_train, n_test, n_features, metric, strategy, dtype, X_train, X_test
        )
