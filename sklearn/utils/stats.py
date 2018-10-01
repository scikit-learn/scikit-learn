import numpy as np

from sklearn.utils.extmath import stable_cumsum


def _weighted_percentile(array, sample_weight, percentile=50):
    """
    Compute the weighted ``percentile`` of ``array`` with ``sample_weight``.
    """
    sorted_idx = np.argsort(array)

    # Find index of median prediction for each sample
    weight_cdf = stable_cumsum(sample_weight[sorted_idx])
    percentile_idx = np.searchsorted(
        weight_cdf, (percentile / 100.) * weight_cdf[-1])
    # in rare cases, percentile_idx equals to len(sorted_idx)
    percentile_idx = np.clip(percentile_idx, 0, len(sorted_idx)-1)
    return array[sorted_idx[percentile_idx]]
