"""
The :mod:`sklearn.manifold` module implements data embedding techniques.
"""

from .locally_linear import locally_linear_embedding, LocallyLinearEmbedding
from .isomap import Isomap
from .minibatch_isomap import MiniBatchIsomap
from .mds import MDS
from .spectral_embedding_ import SpectralEmbedding, spectral_embedding

__all__ = ['locally_linear_embedding', 'LocallyLinearEmbedding', 'Isomap',
           'MiniBatchIsomap', 'MDS', 'SpectralEmbedding', 
           'spectral_embedding']
