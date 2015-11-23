"""Spectral Embedding"""

# Author: Gael Varoquaux <gael.varoquaux@normalesup.org>
#         Wei LI <kuantkid@gmail.com>
# License: BSD 3 clause

from ..utils.deprecation import deprecated

from .laplacian_eigenmap_ import LaplacianEigenmap, laplacian_eigenmap


@deprecated("Function 'spectral_embedding' has been renamed to "
            "'laplacian_eigenmap' and will be removed in release 0.20")
def spectral_embedding(adjacency, n_components=8, eigen_solver=None,
                       random_state=None, eigen_tol=0.0,
                       norm_laplacian=True, drop_first=True):
    return laplacian_eigenmap(adjacency, n_components=n_components,
                              eigen_solver=eigen_solver, drop_first=drop_first,
                              random_state=random_state, eigen_tol=eigen_tol,
                              norm_laplacian=norm_laplacian)
spectral_embedding.__doc__ = laplacian_eigenmap.__doc__


@deprecated("Class 'SpectralEmbedding' has been renamed to 'LaplacianEigenmap'"
            " and will be removed in release 0.20")
class SpectralEmbedding(LaplacianEigenmap):
    pass
