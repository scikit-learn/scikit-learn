"""Data embedding techniques."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from sklearn.manifold._isomap import Isomap
from sklearn.manifold._locally_linear import (
    LocallyLinearEmbedding,
    locally_linear_embedding,
)
from sklearn.manifold._mds import MDS, smacof
from sklearn.manifold._spectral_embedding import SpectralEmbedding, spectral_embedding
from sklearn.manifold._t_sne import TSNE, trustworthiness

__all__ = [
    "MDS",
    "TSNE",
    "Isomap",
    "LocallyLinearEmbedding",
    "SpectralEmbedding",
    "locally_linear_embedding",
    "smacof",
    "spectral_embedding",
    "trustworthiness",
]
