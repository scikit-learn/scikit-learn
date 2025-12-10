"""Utility function to get the number of bootstrap samples."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from numbers import Integral
from warnings import warn


def _get_n_samples_bootstrap(n_samples, max_samples, sample_weight):
    """
    Get the number of samples in a bootstrap sample.

    Parameters
    ----------
    n_samples : int
        Number of samples in the dataset.

    max_samples : None, int or float
        The maximum number of samples to draw.

        - If None, then draw `n_samples` samples.
        - If int, then draw `max_samples` samples.
        - If float, then draw `max_samples * n_samples` unweighted samples or
          `max_samples * sample_weight.sum()` weighted samples.

    sample_weight : array of shape (n_samples,) or None
        Sample weights with frequency semantics when `max_samples` is explicitly
        set to a float or integer value. When keeping the `max_samples=None` default
        value, the equivalence between fitting with integer weighted data points or
        integer repeated data points is no longer guaranteed because the effective
        bootstrap size is no longer guaranteed to be equivalent.

    Returns
    -------
    n_samples_bootstrap : int
        The total number of samples to draw for the bootstrap sample.
    """
    if max_samples is None:
        return n_samples
    elif isinstance(max_samples, Integral):
        return max_samples

    if sample_weight is None:
        weighted_n_samples = n_samples
        weighted_n_samples_msg = f"the number of samples is {weighted_n_samples} "
    else:
        weighted_n_samples = sample_weight.sum()
        weighted_n_samples_msg = (
            f"the total sum of sample weights is {weighted_n_samples} "
        )

    # max_samples Real fractional value relative to weighted_n_samples
    n_samples_bootstrap = max(int(max_samples * weighted_n_samples), 1)
    # Warn when number of bootstrap samples is suspiciously small
    # This heuristic for "suspiciously small" might be adapted if found
    # unsuitable in practice
    if n_samples_bootstrap < max(10, n_samples ** (1 / 3)):
        warn(
            f"Using the fractional value {max_samples=} when {weighted_n_samples_msg}"
            f"results in a low number ({n_samples_bootstrap}) of bootstrap samples. "
            "We recommend passing `max_samples` as an integer instead."
        )
    return n_samples_bootstrap
