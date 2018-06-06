"""
The :mod:`sklearn.manifold` module implements data embedding techniques.
"""

from .locally_linear import locally_linear_embedding, LocallyLinearEmbedding
from .isomap import Isomap
from .mds import MDS, smacof
from .psmds import MDS as PSMDS, pattern_search_mds
from .spectral_embedding_ import SpectralEmbedding, spectral_embedding
from .t_sne import TSNE

__all__ = ['locally_linear_embedding', 'LocallyLinearEmbedding', 'Isomap',
           'MDS', 'smacof', 'PSMDS', 'pattern_search_mds', 'SpectralEmbedding', 'spectral_embedding', "TSNE"]
