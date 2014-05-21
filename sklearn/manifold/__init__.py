"""
The :mod:`sklearn.manifold` module implements data embedding techniques.
"""

from .locally_linear import locally_linear_embedding, LocallyLinearEmbedding
from .isomap import Isomap
from .incremental_isomap import IncrementalIsomap
from .mds import MDS
from .spectral_embedding_ import SpectralEmbedding, spectral_embedding

__all__ = ['locally_linear_embedding', 'LocallyLinearEmbedding', 'Isomap',
           'MDS', 'SpectralEmbedding', 'spectral_embedding']
