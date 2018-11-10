"""Statistical utilities including weighted percentile"""

import numpy as np
from sklearn.utils.extmath import stable_cumsum


def _weighted_percentile(array, sample_weight, percentile=50):
    """Compute the weighted ``percentile`` of ``array``
    with ``sample_weight``.

    This approach follows

               N
        S_N = sum w_k
              k=1

        p_n = 1 / S_N * (x_n - w_n / 2)

        v = v_k + (v_{k + 1} - v_k) * (P - p_k) / (p_{k + 1} - p_k)

    from
    https://en.wikipedia.org/wiki/Percentile#The_weighted_percentile_method.


    Parameters
    ----------
    array : array-like, shape = (n_samples,)
        Array of data on which to calculate the weighted percentile

    sample_weight : array-like, shape = (n_samples,)
        Array of corresponding sample weights with which to calculate
        the weighted percentile

    percentile : int, optional (default: 50)
        Integer value of Pth percentile to compute

    Returns
    -------
    v : float
        Linearly interpolated weighted percentile.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.utils.stats import _weighted_percentile
    >>> weight = np.array([1, 1])
    >>> data = np.array([0, 1])
    >>> _weighted_percentile(data, weight, percentile=0)
    0.0
    >>> _weighted_percentile(data, weight, percentile=50)
    0.5
    >>> _weighted_percentile(data, weight, percentile=90)
    1.0
    """
    if not isinstance(array, np.ndarray):
        array = np.array(array)

    if not isinstance(sample_weight, np.ndarray):
        sample_weight = np.array(sample_weight)

    if (sample_weight < 0).any():
        raise ValueError("sample_weight must contain positive or 0 weights")

    if percentile < 0:
        raise ValueError("percentile must be positive or 0")

    sorted_idx = np.argsort(array)
    sorted_array = array[sorted_idx]

    # if there are no weights, return the min of ``array``
    if sample_weight.sum() == 0:
        return sorted_array[0]

    # Find index of median prediction for each sample
    weight_cdf = stable_cumsum(sample_weight[sorted_idx])
    p_n = 100. / weight_cdf[-1] * (weight_cdf - sample_weight / 2.)
    return np.interp(percentile, p_n, sorted_array)
