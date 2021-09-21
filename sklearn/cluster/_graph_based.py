"""Graph-Based clustering using connected components and spanning trees."""

# Authors: Dani El-Ayyass <dayyass@yandex.ru>
# License: BSD 3 clause


from typing import Callable, Optional, Tuple, Union

import numpy as np

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
