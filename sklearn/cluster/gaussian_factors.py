""" Algorithms for clustering : Meanshift,  Affinity propagation and spectral
clustering.

"""
# Author: x0l <x0l@github.com>
# Licence: BSD 3 clause

import numpy as np

from ..base import BaseEstimator, ClusterMixin
from ..utils import as_float_array, check_array
from ..utils.validation import check_is_fitted

from . import _gaussian_factors

"""
# example
import numpy as np
from sklearn.cluster import gaussian_factors
from sklearn.datasets import make_blobs

X, y = make_blobs(n_samples=100, n_features=1000)
C = np.corrcoef(X)

Z, labels = gaussian_factors(C)
assert(labels.max()+1 == 3)
"""

def gaussian_factors(X, damping=1e-4):
    """
    Performs Gaussian factors likelihood clustering of data using an
    agglomerative hierarchical algorithm.

    Parameters
    ----------
    X : array, shape (n_samples, n_samples)
        Correlation matrix of samples

    damping : float, optional, default: 1e-4
        Damping factor between 0 and 1.

    Returns
    -------
    Z : array, shape (n_samples - 1, 4)
        The hierarchical clustering encoded as a linkage matrix. See
        `scipy.cluster.hierarchy.linkage` for a description. Here, distances
        are replaced by log-likelihoods (which may be negative) hence a
        preprocessing step is required before using functions like
        `scipy.cluster.hierarchy.dendrogram`.

    labels : array, shape (n_samples,)
        Cluster labels for each point.
    """
    return _gaussian_factors.gaussian_factors(X, damping)


###############################################################################


class GaussianFactors(BaseEstimator, ClusterMixin):
    """Perform Affinity Propagation Clustering of data.

    Parameters
    ----------
    damping : float, optional, default: 1e-4
        Damping factor between 0 and 1.

    max_iter : int, optional, default: 200
        Maximum number of iterations.

    verbose : boolean, optional, default: False
        Whether to be verbose.


    Attributes
    ----------
    labels_ : array, shape (n_samples,)
        Labels of each point

    n_iter_ : int
        Number of iterations taken to converge.

    Notes
    -----
    See examples/cluster/plot_affinity_propagation.py for an example.

    The algorithmic complexity of affinity propagation is quadratic
    in the number of points.

    References
    ----------

    Lorenzo Giada and Matteo Marsili, "Data clustering and noise undressing
    for correlation matrices", Phys. Rev. E 63, 061101 (2001)
    """

    def __init__(self, damping=1e-4, max_iter=200, verbose=False):

        self.damping = damping
        self.max_iter = max_iter
        self.verbose = verbose


    def fit(self, X):
        # TODO save what's needed to perform prediction...
        pass

    def predict(self, X):
        # TODO use max likelihood
        pass
