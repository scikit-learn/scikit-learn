"""
The :mod:`sklearn.manifold` module implements data embedding techniques.
"""

from ._graph_spectral_embedding import GraphSpectralEmbedding
from ._locally_linear import locally_linear_embedding, LocallyLinearEmbedding
from ._isomap import Isomap
from ._mds import MDS, smacof
from ._spectral_embedding import SpectralEmbedding, spectral_embedding
from ._t_sne import TSNE, trustworthiness

__all__ = [
    "GraphSpectralEmbedding",
    "locally_linear_embedding",
    "LocallyLinearEmbedding",
    "Isomap",
    "MDS",
    "smacof",
    "SpectralEmbedding",
    "spectral_embedding",
    "TSNE",
    "trustworthiness",
]
