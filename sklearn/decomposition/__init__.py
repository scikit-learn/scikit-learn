"""
The :mod:`sklearn.decomposition` module includes matrix decomposition
algorithms, including among others PCA, NMF or ICA. Most of the algorithms of
this module can be regarded as dimensionality reduction techniques.
"""
from ..externals import _lazy_loader

lazy__getattr__, lazy__dir__, lazy__all__ = _lazy_loader.attach(
    __name__,
    submod_attrs={
        "_dict_learning": [
            "DictionaryLearning",
            "MiniBatchDictionaryLearning",
            "SparseCoder",
            "dict_learning",
            "dict_learning_online",
            "sparse_encode",
        ],
        "_factor_analysis": ["FactorAnalysis"],
        "_fastica": ["FastICA", "fastica"],
        "_incremental_pca": ["IncrementalPCA"],
        "_kernel_pca": ["KernelPCA"],
        "_lda": ["LatentDirichletAllocation"],
        "_nmf": ["MiniBatchNMF", "NMF", "non_negative_factorization"],
        "_pca": ["PCA"],
        "_sparse_pca": ["MiniBatchSparsePCA", "SparsePCA"],
        "_truncated_svd": ["TruncatedSVD"],
    },
)

# randomized_svd is importable but is defined in utils.extmath
__all__ = lazy__all__ + ["randomized_svd"]


def __dir__():
    return __all__


def __getattr__(name):
    if name == "randomized_svd":
        from ..utils.extmath import randomized_svd

        return randomized_svd
    return lazy__getattr__(name)
