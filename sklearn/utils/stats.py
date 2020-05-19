import numpy as np

from .extmath import stable_cumsum
from sklearn.utils.fixes import _take_along_axis


def _weighted_percentile(array, sample_weight, percentile=50):
    """
    Compute the weighted ``percentile`` of ``array`` with ``sample_weight``.
    If ``array`` is 2D, compute weighted ``percentile`` along axis=0.
    """
    n_dim = array.ndim
    if n_dim == 0:
        return array[()]
    if array.ndim == 1:
        array = array.reshape((-1, 1))
    if array.shape != sample_weight.shape:
        sample_weight = sample_weight.reshape(array.shape)
    sorted_idx = np.argsort(array, axis=0)
    sorted_weights = _take_along_axis(sample_weight, sorted_idx, axis=0)

    # Find index of median prediction for each sample
    weight_cdf = stable_cumsum(sorted_weights, axis=0)
    adjusted_percentile = (percentile / 100.) * weight_cdf[-1]
    percentile_idx = [np.searchsorted(weight_cdf[:, i], adjusted_percentile[i])
                      for i in range(weight_cdf.shape[1])]
    percentile_idx = np.array(percentile_idx)
    # In rare cases, percentile_idx equals to sorted_idx.shape[0]
    max_idx = sorted_idx.shape[0] - 1
    percentile_idx = np.apply_along_axis(lambda x: np.clip(x, 0, max_idx),
                                         axis=0, arr=percentile_idx)

    col_index = np.arange(array.shape[1])
    percentile_in_sorted = sorted_idx[percentile_idx, col_index]
    percentile = array[percentile_in_sorted, col_index]
    return percentile[0] if ndim == 1 else percentile
