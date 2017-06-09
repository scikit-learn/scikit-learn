# Remove in version 0.21

from scipy.sparse.csgraph import connected_components as \
     scipy_connected_components

from sklearn.utils.deprecation import deprecated


@deprecated("sklearn.utils.sparsetools.connected_components was deprecated in "
            "version 0.19 and will be removed in 0.21. Use "
            "scipy.sparse.csgraph.connected_components instead.")
def connected_components(*args, **kwargs):
    return scipy_connected_components(*args, **kwargs)
