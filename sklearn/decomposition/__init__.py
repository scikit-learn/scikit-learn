"""Matrix decomposition algorithms.

These include PCA, NMF, ICA, and more. Most of the algorithms of this module can be
regarded as dimensionality reduction techniques.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from sklearn.decomposition._dict_learning import (
    DictionaryLearning,
    MiniBatchDictionaryLearning,
    SparseCoder,
    dict_learning,
    dict_learning_online,
    sparse_encode,
)
from sklearn.decomposition._factor_analysis import FactorAnalysis
from sklearn.decomposition._fastica import FastICA, fastica
from sklearn.decomposition._incremental_pca import IncrementalPCA
from sklearn.decomposition._kernel_pca import KernelPCA
from sklearn.decomposition._lda import LatentDirichletAllocation
from sklearn.decomposition._nmf import NMF, MiniBatchNMF, non_negative_factorization
from sklearn.decomposition._pca import PCA
from sklearn.decomposition._sparse_pca import MiniBatchSparsePCA, SparsePCA
from sklearn.decomposition._truncated_svd import TruncatedSVD
from sklearn.utils.extmath import randomized_svd

__all__ = [
    "NMF",
    "PCA",
    "DictionaryLearning",
    "FactorAnalysis",
    "FastICA",
    "IncrementalPCA",
    "KernelPCA",
    "LatentDirichletAllocation",
    "MiniBatchDictionaryLearning",
    "MiniBatchNMF",
    "MiniBatchSparsePCA",
    "SparseCoder",
    "SparsePCA",
    "TruncatedSVD",
    "dict_learning",
    "dict_learning_online",
    "fastica",
    "non_negative_factorization",
    "randomized_svd",
    "sparse_encode",
]
