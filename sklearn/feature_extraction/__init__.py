"""
The :mod:`sklearn.feature_extraction` module deals with feature extraction
from raw data. It currently includes methods to extract features from text and
images.
"""
from ..externals import _lazy_loader

__getattr__, __dir__, __all__ = _lazy_loader.attach(
    __name__,
    ["image", "text"],
    submod_attrs={
        "_dict_vectorizer": ["DictVectorizer"],
        "_hash": ["FeatureHasher"],
        "image": ["img_to_graph", "grid_to_graph"],
    },
)
