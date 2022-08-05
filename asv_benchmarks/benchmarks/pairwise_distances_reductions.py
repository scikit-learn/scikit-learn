import numpy as np
from scipy.spatial.distance import cdist

from .common import Benchmark

from sklearn.metrics._pairwise_distances_reduction import (
    PairwiseDistancesArgKmin,
    PairwiseDistancesRadiusNeighborhood,
)

# To run benchmarks defined this file, between for instance your <current_branch>
# and its fork point with upstream/main, run:
#
#             asv continuous -b PairwiseDistancesReductions \
#             -e `git merge-base --fork-point upstream/main <your_branch>` <your_branch>


class PairwiseDistancesReductionsBenchmark(Benchmark):

    param_names = ["n_train", "n_test", "n_features", "metric", "strategy"]
    params = [
        [1000, 10_000, int(1e7)],
        [1000, 10_000, 100_000],
        [100],
        ["euclidean", "manhattan"],
        ["auto", "parallel_on_X", "parallel_on_Y"],
    ]

    def setup(self, n_train, n_test, n_features, metric, strategy):
        rng = np.random.RandomState(0)
        self.X_train = rng.rand(n_train, n_features)
        self.X_test = rng.rand(n_test, n_features)
        self.y_train = rng.randint(low=-1, high=1, size=(n_train,))
        self.metric = metric
        self.strategy = strategy

        self.k = 10

        # Motive: radius has to be scaled with the number of feature
        # Defining it as the 0.001-quantile allows to have in expectation
        # a constant amount of neighbors, regardless of the value of n_features.
        self.radius = np.quantile(
            a=cdist(self.X_train[:10000], self.X_test[:10]).ravel(), q=0.001
        )

    def time_PairwiseDistancesArgKmin(
        self, n_train, n_test, n_features, metric, strategy
    ):

        PairwiseDistancesArgKmin.compute(
            X=self.X_test,
            Y=self.X_train,
            k=10,
            metric=self.metric,
            return_distance=True,
            strategy=self.strategy,
        )

    def peakmem_PairwiseDistancesArgKmin(
        self, n_train, n_test, n_features, metric, strategy
    ):
        self.time_PairwiseDistancesArgKmin(
            n_train, n_test, n_features, metric, strategy
        )

    def time_PairwiseDistancesRadiusNeighborhood(
        self, n_train, n_test, n_features, metric, strategy
    ):

        PairwiseDistancesRadiusNeighborhood.compute(
            X=self.X_test,
            Y=self.X_train,
            radius=self.radius,
            metric=self.metric,
            return_distance=True,
            strategy=self.strategy,
        )

    def peakmem_PairwiseDistancesRadiusNeighborhood(
        self, n_train, n_test, n_features, metric, strategy
    ):
        self.time_PairwiseDistancesRadiusNeighborhood(
            n_train, n_test, n_features, metric, strategy
        )
