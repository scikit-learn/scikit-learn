"""
    sphinx.search.ro
    ~~~~~~~~~~~~~~~~

    Romanian search language: includes the JS Romanian stemmer.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

from typing import Dict, Set

import snowballstemmer

from sphinx.search import SearchLanguage


class SearchRomanian(SearchLanguage):
    lang = 'ro'
    language_name = 'Romanian'
    js_stemmer_rawcode = 'romanian-stemmer.js'
    stopwords: Set[str] = set()

    def init(self, options: Dict) -> None:
        self.stemmer = snowballstemmer.stemmer('romanian')

    def stem(self, word: str) -> str:
        return self.stemmer.stemWord(word.lower())
