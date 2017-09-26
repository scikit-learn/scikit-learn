import numpy as np
from scipy.stats import rankdata as scipy_rankdata

from sklearn.utils.extmath import stable_cumsum
from sklearn.utils.deprecation import deprecated


# Remove in sklearn 0.21
@deprecated("sklearn.utils.stats.rankdata was deprecated in version 0.19 and "
            "will be removed in 0.21. Use scipy.stats.rankdata instead.")
def rankdata(*args, **kwargs):
    return scipy_rankdata(*args, **kwargs)


def _weighted_percentile(array, sample_weight, percentile=50):
    """
    Compute the weighted ``percentile`` of ``array`` with ``sample_weight``.
    """
    sorted_idx = np.argsort(array)

    # Find index of median prediction for each sample
    weight_cdf = stable_cumsum(sample_weight[sorted_idx])
    percentile_idx = np.searchsorted(
        weight_cdf, (percentile / 100.) * weight_cdf[-1])
    return array[sorted_idx[percentile_idx]]
