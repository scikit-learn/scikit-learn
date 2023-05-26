"""
The :mod:`sklearn.manifold` module implements data embedding techniques.
"""
from ..externals import _lazy_loader

__getattr__, __dir__, __all__ = _lazy_loader.attach(
    __name__,
    submod_attrs={
        "_isomap": ["Isomap"],
        "_locally_linear": ["locally_linear_embedding", "LocallyLinearEmbedding"],
        "_mds": ["MDS", "smacof"],
        "_spectral_embedding": ["SpectralEmbedding", "spectral_embedding"],
        "_t_sne": ["TSNE", "trustworthiness"],
    },
)
