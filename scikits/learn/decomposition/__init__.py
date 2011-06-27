"""
Matrix decomposition algorithms
"""

from .nmf import NMF, ProjectedGradientNMF
from .pca import PCA, RandomizedPCA, ProbabilisticPCA, KernelPCA
from .sparse_pca import SparsePCA, dict_learning
from .fastica_ import FastICA, fastica
