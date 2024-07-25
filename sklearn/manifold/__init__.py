"""Data embedding techniques."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from ._isomap import Isomap
from ._locally_linear import LocallyLinearEmbedding, locally_linear_embedding
from ._mds import MDS, smacof
from ._spectral_embedding import SpectralEmbedding, spectral_embedding
from ._t_sne import TSNE, trustworthiness

__all__ = [
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
