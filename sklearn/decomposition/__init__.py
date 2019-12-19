"""
The :mod:`sklearn.decomposition` module includes matrix decomposition
algorithms, including among others PCA, NMF or ICA. Most of the algorithms of
this module can be regarded as dimensionality reduction techniques.
"""

# TODO: remove me in 0.24 (as well as the noqa markers) and
# import the dict_learning func directly from the ._dict_learning
# module instead.
# Pre-cache the import of the deprecated module so that import
# sklearn.decomposition.dict_learning returns the function as in
# 0.21, instead of the module.
# https://github.com/scikit-learn/scikit-learn/issues/15842
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=FutureWarning)
    from .dict_learning import dict_learning


from ._nmf import NMF, non_negative_factorization  # noqa
from ._pca import PCA  # noqa
from ._incremental_pca import IncrementalPCA  # noqa
from ._kernel_pca import KernelPCA  # noqa
from ._sparse_pca import SparsePCA, MiniBatchSparsePCA  # noqa
from ._truncated_svd import TruncatedSVD  # noqa
from ._fastica import FastICA, fastica  # noqa
from ._dict_learning import (dict_learning_online,
                             sparse_encode, DictionaryLearning,
                             MiniBatchDictionaryLearning, SparseCoder)  # noqa
from ._factor_analysis import FactorAnalysis  # noqa
from ..utils.extmath import randomized_svd  # noqa
from ._lda import LatentDirichletAllocation  # noqa


__all__ = ['DictionaryLearning',
           'FastICA',
           'IncrementalPCA',
           'KernelPCA',
           'MiniBatchDictionaryLearning',
           'MiniBatchSparsePCA',
           'NMF',
           'PCA',
           'SparseCoder',
           'SparsePCA',
           'dict_learning',
           'dict_learning_online',
           'fastica',
           'non_negative_factorization',
           'randomized_svd',
           'sparse_encode',
           'FactorAnalysis',
           'TruncatedSVD',
           'LatentDirichletAllocation']
