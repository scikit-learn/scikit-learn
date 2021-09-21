"""Graph-Based clustering using connected components and spanning trees."""

# Authors: Dani El-Ayyass <dayyass@yandex.ru>
# License: BSD 3 clause


from typing import Callable, Optional, Tuple, Union

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from ..base import BaseEstimator, ClusterMixin
from ..metrics import pairwise_distances


def _pairwise_distances(
    X: np.ndarray,
    metric: Union[str, Callable] = "euclidean",
    n_jobs: Optional[int] = None,
) -> np.ndarray:
    """TODO"""

    distances = pairwise_distances(X=X, metric=metric, n_jobs=n_jobs)

    return distances


def distances_to_adjacency_matrix(
    distances: np.ndarray,
    threshold: float,
) -> np.ndarray:
    """TODO"""

    N = distances.shape[0]

    adjacency_matrix = (distances < threshold).astype(int) - np.eye(N, dtype=int)

    return adjacency_matrix


def span_tree_top_n_weights_idx(
    span_tree: np.ndarray,
    n: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """TODO"""

    span_tree_shape = span_tree.shape
    N = span_tree_shape[0]

    if n == 0:
        array_1_to_N = np.array(range(N))
        unravel_idx = (array_1_to_N, array_1_to_N)

    else:
        ravel_span_tree_top_n_weights_idx = np.argpartition(
            a=span_tree.ravel(),
            kth=-n,
        )[-n:]

        unravel_idx = np.unravel_index(
            indices=ravel_span_tree_top_n_weights_idx,
            shape=span_tree_shape,
        )

    return unravel_idx


class ConnectedComponentsClustering(ClusterMixin, BaseEstimator):

    """TODO"""

    def __init__(
        self,
        threshold: float,
        metric: Union[str, Callable] = "euclidean",
        n_jobs: Optional[int] = None,
    ) -> None:
        """TODO"""

        self.threshold = threshold
        self.metric = metric
        self.n_jobs = n_jobs

    def fit(
        self,
        X: np.ndarray,
        y: np.ndarray = None,
    ):
        """TODO"""

        X = self._validate_data(X, accept_sparse="csr")

        distances = _pairwise_distances(
            X=X,
            metric=self.metric,
            n_jobs=self.n_jobs,
        )

        adjacency_matrix = distances_to_adjacency_matrix(
            distances=distances,
            threshold=self.threshold,
        )

        graph = csr_matrix(adjacency_matrix)

        _, labels = connected_components(
            csgraph=graph,
            directed=True,
            return_labels=True,
        )

        self.labels_ = labels

        return self

    def fit_predict(
        self,
        X: np.ndarray,
        y: np.ndarray = None,
    ):
        """TODO"""

        self.fit(X)
        return self.labels_
