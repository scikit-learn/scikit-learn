"""
The :mod:`sklearn.manifold` module implements data embedding techniques.
"""

from .locally_linear import locally_linear_embedding, LocallyLinearEmbedding
from .isomap import Isomap
from .mds import MDS, smacof
from .spectral_embedding_ import SpectralEmbedding, spectral_embedding
from .t_sne import TSNE

__all__ = ['locally_linear_embedding', 'LocallyLinearEmbedding', 'Isomap',
           'MDS', 'smacof', 'SpectralEmbedding', 'spectral_embedding', "TSNE"]
