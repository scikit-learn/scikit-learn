"""
The :mod:`sklearn.feature_extraction` module deals with feature extraction
from raw data. It currently includes methods to extract features from text and
images.
"""

from .dict_vectorizer import DictVectorizer
from .hashing import FeatureHasher
from .image import img_to_graph, grid_to_graph
from . import text

from . import tokenize_utils
from . import syllables_en
from . import readability

__all__ = ['DictVectorizer', 'image', 'img_to_graph', 'grid_to_graph', 'text',
           'FeatureHasher',
           'tokenize_utils', 'syllables_en', 'readability']
