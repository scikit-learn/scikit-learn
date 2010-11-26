"""Utilities to preprocess text content and vectorize it

The vectorizers are able to output both dense and sparse representations based
on the implementation used.
"""

from .dense import CharNGramAnalyzer
from .dense import CountVectorizer
from .dense import DEFAULT_ANALYZER
from .dense import DEFAULT_PREPROCESSOR
from .dense import ENGLISH_STOP_WORDS
from .dense import RomanPreprocessor
from .dense import TfidfTransformer
from .dense import Vectorizer
from .dense import WordNGramAnalyzer
from .dense import strip_accents
from .dense import strip_tags
from .dense import to_ascii
