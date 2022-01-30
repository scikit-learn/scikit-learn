"""
    sphinx.util.stemmer
    ~~~~~~~~~~~~~~~~~~~

    Word stemming utilities for Sphinx.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

from sphinx.util.stemmer.porter import PorterStemmer

try:
    from Stemmer import Stemmer as _PyStemmer
    PYSTEMMER = True
except ImportError:
    PYSTEMMER = False


class BaseStemmer:
    def stem(self, word: str) -> str:
        raise NotImplementedError()


class PyStemmer(BaseStemmer):
    def __init__(self) -> None:
        self.stemmer = _PyStemmer('porter')

    def stem(self, word: str) -> str:
        return self.stemmer.stemWord(word)


class StandardStemmer(PorterStemmer, BaseStemmer):
    """All those porter stemmer implementations look hideous;
    make at least the stem method nicer.
    """
    def stem(self, word: str) -> str:  # type: ignore
        return super().stem(word, 0, len(word) - 1)


def get_stemmer() -> BaseStemmer:
    if PYSTEMMER:
        return PyStemmer()
    else:
        return StandardStemmer()
