"""
The :mod:`sklearn.decomposition` module includes matrix decomposition
algorithms, including among others PCA, NMF or ICA. Most of the algorithms of
this module can be regarded as dimensionality reduction techniques.
"""

from .nmf import NMF, ProjectedGradientNMF
from .pca import PCA, RandomizedPCA, ProbabilisticPCA
from .kernel_pca import KernelPCA
from .sparse_pca import SparsePCA, MiniBatchSparsePCA
from .fastica_ import FastICA, fastica
from .dict_learning import dict_learning, dict_learning_online, sparse_encode,\
                           DictionaryLearning, MiniBatchDictionaryLearning,\
                           SparseCoder

__all__ = ['DictionaryLearning',
           'FastICA',
           'KernelPCA',
           'MiniBatchDictionaryLearning',
           'MiniBatchSparsePCA',
           'NMF',
           'PCA',
           'ProbabilisticPCA',
           'ProjectedGradientNMF',
           'RandomizedPCA',
           'SparseCoder',
           'SparsePCA',
           'dict_learning',
           'dict_learning_online',
           'fastica',
           'sparse_encode']
