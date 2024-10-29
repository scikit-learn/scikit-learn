"""Matrix decomposition algorithms.

These include PCA, NMF, ICA, and more. Most of the algorithms of this module can be
regarded as dimensionality reduction techniques.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from ..utils.extmath import randomized_svd
from ._dict_learning import (
    DictionaryLearning,
    MiniBatchDictionaryLearning,
    SparseCoder,
    dict_learning,
    dict_learning_online,
    sparse_encode,
)
from ._factor_analysis import FactorAnalysis
from ._fastica import FastICA, fastica
from ._incremental_pca import IncrementalPCA
from ._kernel_pca import KernelPCA
from ._lda import LatentDirichletAllocation
from ._nmf import (
    NMF,
    MiniBatchNMF,
    non_negative_factorization,
)
from ._pca import PCA
from ._sparse_pca import MiniBatchSparsePCA, SparsePCA
from ._truncated_svd import TruncatedSVD

__all__ = [
    "DictionaryLearning",
    "FastICA",
    "IncrementalPCA",
    "KernelPCA",
    "MiniBatchDictionaryLearning",
    "MiniBatchNMF",
    "MiniBatchSparsePCA",
    "NMF",
    "PCA",
    "SparseCoder",
    "SparsePCA",
    "dict_learning",
    "dict_learning_online",
    "fastica",
    "non_negative_factorization",
    "randomized_svd",
    "sparse_encode",
    "FactorAnalysis",
    "TruncatedSVD",
    "LatentDirichletAllocation",
]
