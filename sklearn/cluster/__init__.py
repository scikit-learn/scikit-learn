"""
The :mod:`sklearn.cluster` module gathers popular unsupervised clustering
algorithms.
"""

from .spectral import spectral_clustering, SpectralClustering
from .mean_shift_ import mean_shift, MeanShift, estimate_bandwidth, \
    get_bin_seeds
from .affinity_propagation_ import affinity_propagation, AffinityPropagation
from .hierarchical import ward_tree, Ward, WardAgglomeration
from .k_means_ import k_means, KMeans, MiniBatchKMeans
from .dbscan_ import dbscan, DBSCAN
from .bicluster import SpectralBiclustering, SpectralCoclustering
from ..utils import deprecated


# backward compatibility
@deprecated("to be removed in 0.15;"
            " use sklearn.manifold.spectral_embedding instead")
def spectral_embedding(*args, **kwargs):
    """Deprecated, use ``sklearn.manifold.spectral_embedding`` instead"""
    from ..manifold.spectral_embedding_ import spectral_embedding
    return spectral_embedding(*args, **kwargs)


__all__ = ['AffinityPropagation',
           'DBSCAN',
           'KMeans',
           'MeanShift',
           'MiniBatchKMeans',
           'SpectralClustering',
           'Ward',
           'WardAgglomeration',
           'affinity_propagation',
           'dbscan',
           'estimate_bandwidth',
           'get_bin_seeds',
           'k_means',
           'mean_shift',
           'spectral_clustering',
           'spectral_embedding',
           'ward_tree',
           'SpectralBiclustering',
           'SpectralCoclustering']
