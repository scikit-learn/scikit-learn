"""
Matrix decomposition algorithms
"""

from .nmf import NMF, ProjectedGradientNMF
from .pca import PCA, RandomizedPCA, ProbabilisticPCA
from .kernel_pca import KernelPCA
from .sparse_pca import SparsePCA, MiniBatchSparsePCA
from .fastica_ import FastICA, fastica
from .dict_learning import dict_learning, dict_learning_online, \
                           DictionaryLearning, MiniBatchDictionaryLearning, \
                           sparse_encode, sparse_encode_parallel
