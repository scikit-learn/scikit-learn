"""
The :mod:`sklearn.manifold` module implements data embedding techniques.
"""

from .locally_linear import locally_linear_embedding, LocallyLinearEmbedding
from .isomap import Isomap
from .mds import MDS
from .spectral_embedding_ import SpectralEmbedding, spectral_embedding
from .t_sne import TSNE
from .laplacian_eigenmap_ import LaplacianEigenmap, laplacian_eigenmap

__all__ = ['locally_linear_embedding', 'LocallyLinearEmbedding', 'Isomap',
           'MDS', "TSNE", "LaplacianEigenmap", "laplacian_eigenmap"]
