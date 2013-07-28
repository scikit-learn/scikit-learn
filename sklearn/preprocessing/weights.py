import numpy as np

from sklearn.utils import safe_asarray


def balance_weights(y):
    """Compute sample weights such that the class distribution of y becomes
       balanced.

    Parameters
    ----------
    y : array-like
        Labels for the samples.

    Returns
    -------
    weights : array-like
        The sample weights.
    """
    y = safe_asarray(y)
    y = np.searchsorted(np.unique(y), y)
    bins = np.bincount(y)

    weights = 1. / bins.take(y)
    weights *= bins.min()

    return weights
