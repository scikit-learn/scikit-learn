import numpy as np

from .extmath import stable_cumsum
from sklearn.utils.fixes import _take_along_axis


def _weighted_percentile(array, sample_weight, percentile=50):
    """
    Compute the weighted ``percentile`` of ``array`` with ``sample_weight``.
    If ``array`` is 2D, compute weighted ``percentile`` along axis=0.
    """
    if sample_weight is None:
        sample_weight = np.ones_like(array)

    n_dim = array.ndim
    if n_dim == 0:
        return array[()]
    # Remove trailing single-dimensional entries from the array shape
    if n_dim == 2 and array.shape[1] == 1:
        array.reshape((array.shape[0],))
    sorted_idx = np.argsort(array, axis=0)
    sorted_weights = _take_along_axis(sample_weight, sorted_idx, axis=0)

    # Find index of median prediction for each sample
    weight_cdf = stable_cumsum(sorted_weights, axis=0)
    adjusted_percentile = (percentile / 100.) * weight_cdf[-1]
    if n_dim == 1:
        percentile_idx = [np.searchsorted(weight_cdf, adjusted_percentile)]
    elif n_dim == 2:
        percentile_idx = [np.searchsorted(weight_cdf[:, i],
                                          adjusted_percentile[i])
                          for i in range(weight_cdf.shape[1])]
    percentile_idx = np.array(percentile_idx)
    # In rare cases, percentile_idx equals to sorted_idx.shape[0]
    max_idx = sorted_idx.shape[0] - 1
    percentile_idx = np.apply_along_axis(lambda x: np.clip(x, 0, max_idx),
                                         axis=0, arr=percentile_idx)
    if n_dim == 1:
        percentile = array[sorted_idx[percentile_idx]][0]
    elif n_dim == 2:
        col_index = np.arange(array.shape[1])
        fancy_index_percentile = tuple([percentile_idx, col_index])
        percentile_in_sorted = sorted_idx[fancy_index_percentile]
        fancy_index_sorted = tuple([percentile_in_sorted, col_index])
        percentile = array[fancy_index_sorted]
    return percentile
