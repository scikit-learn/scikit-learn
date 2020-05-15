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
    array = np.squeeze(array)
    sample_weight = np.squeeze(sample_weight)

    n_dim = array.ndim
    if n_dim == 0:
        return array[()]
    sorted_idx = np.argsort(array, axis=0)
    sorted_weights = _take_along_axis(sample_weight, sorted_idx, axis=0)
    sorted_array = _take_along_axis(array, sorted_idx, axis=0)

    # Find index of median prediction for each sample
    weight_cdf = stable_cumsum(sorted_weights, axis=0)
    max_weight_cdf = np.take(weight_cdf, -1, axis=0)
    adjusted_percentile = (percentile / 100.) * max_weight_cdf
    if n_dim == 1:
        percentile_idx = [np.searchsorted(weight_cdf, adjusted_percentile)]
    elif n_dim == 2:
        percentile_idx = [np.searchsorted(weight_cdf[:, i],
                                          adjusted_percentile[i])
                          for i in range(weight_cdf.shape[1])]
    percentile_idx = np.array(percentile_idx)
    # in rare cases, percentile_idx equals to sorted_idx.shape[0]
    max_idx = sorted_idx.shape[0] - 1
    percentile_idx = np.apply_along_axis(lambda x: np.clip(x, 0, max_idx),
                                         axis=0, arr=percentile_idx)
    if n_dim == 1:
        percentile = sorted_array[percentile_idx][0]
    elif n_dim == 2:
        n_col = sorted_array.shape[1]
        percentile = sorted_array[percentile_idx, np.arange(n_col)]
    return percentile
