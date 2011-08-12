"""
Matrix decomposition algorithms
"""

from .nmf import NMF, ProjectedGradientNMF
from .pca import PCA, RandomizedPCA, ProbabilisticPCA
from .kernel_pca import KernelPCA
from .sparse_pca import SparsePCA, MiniBatchSparsePCA, dict_learning, \
                        dict_learning_online
from .fastica_ import FastICA, fastica
from .dict_learning import DictionaryLearning, DictionaryLearningOnline
