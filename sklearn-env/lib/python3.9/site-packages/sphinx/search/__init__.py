"""
    sphinx.search
    ~~~~~~~~~~~~~

    Create a full-text search index for offline search.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""
import html
import pickle
import re
from importlib import import_module
from os import path
from typing import IO, Any, Dict, Iterable, List, Optional, Set, Tuple, Type

from docutils import nodes
from docutils.nodes import Element, Node

from sphinx import addnodes, package_dir
from sphinx.environment import BuildEnvironment
from sphinx.search.jssplitter import splitter_code
from sphinx.util import jsdump


class SearchLanguage:
    """
    This class is the base class for search natural language preprocessors.  If
    you want to add support for a new language, you should override the methods
    of this class.

    You should override `lang` class property too (e.g. 'en', 'fr' and so on).

    .. attribute:: stopwords

       This is a set of stop words of the target language.  Default `stopwords`
       is empty.  This word is used for building index and embedded in JS.

    .. attribute:: js_splitter_code

       Return splitter function of JavaScript version.  The function should be
       named as ``splitQuery``.  And it should take a string and return list of
       strings.

       .. versionadded:: 3.0

    .. attribute:: js_stemmer_code

       Return stemmer class of JavaScript version.  This class' name should be
       ``Stemmer`` and this class must have ``stemWord`` method.  This string is
       embedded as-is in searchtools.js.

       This class is used to preprocess search word which Sphinx HTML readers
       type, before searching index. Default implementation does nothing.
    """
    lang: str = None
    language_name: str = None
    stopwords: Set[str] = set()
    js_splitter_code: str = None
    js_stemmer_rawcode: str = None
    js_stemmer_code = """
/**
 * Dummy stemmer for languages without stemming rules.
 */
var Stemmer = function() {
  this.stemWord = function(w) {
    return w;
  }
}
"""

    _word_re = re.compile(r'(?u)\w+')

    def __init__(self, options: Dict) -> None:
        self.options = options
        self.init(options)

    def init(self, options: Dict) -> None:
        """
        Initialize the class with the options the user has given.
        """

    def split(self, input: str) -> List[str]:
        """
        This method splits a sentence into words.  Default splitter splits input
        at white spaces, which should be enough for most languages except CJK
        languages.
        """
        return self._word_re.findall(input)

    def stem(self, word: str) -> str:
        """
        This method implements stemming algorithm of the Python version.

        Default implementation does nothing.  You should implement this if the
        language has any stemming rules.

        This class is used to preprocess search words before registering them in
        the search index.  The stemming of the Python version and the JS version
        (given in the js_stemmer_code attribute) must be compatible.
        """
        return word

    def word_filter(self, word: str) -> bool:
        """
        Return true if the target word should be registered in the search index.
        This method is called after stemming.
        """
        return (
            len(word) == 0 or not (
                ((len(word) < 3) and (12353 < ord(word[0]) < 12436)) or
                (ord(word[0]) < 256 and (
                    word in self.stopwords
                ))))


# SearchEnglish imported after SearchLanguage is defined due to circular import
from sphinx.search.en import SearchEnglish


def parse_stop_word(source: str) -> Set[str]:
    """
    Parse snowball style word list like this:

    * http://snowball.tartarus.org/algorithms/finnish/stop.txt
    """
    result: Set[str] = set()
    for line in source.splitlines():
        line = line.split('|')[0]  # remove comment
        result.update(line.split())
    return result


# maps language name to module.class or directly a class
languages: Dict[str, Any] = {
    'da': 'sphinx.search.da.SearchDanish',
    'de': 'sphinx.search.de.SearchGerman',
    'en': SearchEnglish,
    'es': 'sphinx.search.es.SearchSpanish',
    'fi': 'sphinx.search.fi.SearchFinnish',
    'fr': 'sphinx.search.fr.SearchFrench',
    'hu': 'sphinx.search.hu.SearchHungarian',
    'it': 'sphinx.search.it.SearchItalian',
    'ja': 'sphinx.search.ja.SearchJapanese',
    'nl': 'sphinx.search.nl.SearchDutch',
    'no': 'sphinx.search.no.SearchNorwegian',
    'pt': 'sphinx.search.pt.SearchPortuguese',
    'ro': 'sphinx.search.ro.SearchRomanian',
    'ru': 'sphinx.search.ru.SearchRussian',
    'sv': 'sphinx.search.sv.SearchSwedish',
    'tr': 'sphinx.search.tr.SearchTurkish',
    'zh': 'sphinx.search.zh.SearchChinese',
}


class _JavaScriptIndex:
    """
    The search index as JavaScript file that calls a function
    on the documentation search object to register the index.
    """

    PREFIX = 'Search.setIndex('
    SUFFIX = ')'

    def dumps(self, data: Any) -> str:
        return self.PREFIX + jsdump.dumps(data) + self.SUFFIX

    def loads(self, s: str) -> Any:
        data = s[len(self.PREFIX):-len(self.SUFFIX)]
        if not data or not s.startswith(self.PREFIX) or not \
           s.endswith(self.SUFFIX):
            raise ValueError('invalid data')
        return jsdump.loads(data)

    def dump(self, data: Any, f: IO) -> None:
        f.write(self.dumps(data))

    def load(self, f: IO) -> Any:
        return self.loads(f.read())


js_index = _JavaScriptIndex()


class WordCollector(nodes.NodeVisitor):
    """
    A special visitor that collects words for the `IndexBuilder`.
    """

    def __init__(self, document: nodes.document, lang: SearchLanguage) -> None:
        super().__init__(document)
        self.found_words: List[str] = []
        self.found_title_words: List[str] = []
        self.lang = lang

    def is_meta_keywords(self, node: Element) -> bool:
        if (isinstance(node, (addnodes.meta, addnodes.docutils_meta)) and
                node.get('name') == 'keywords'):
            meta_lang = node.get('lang')
            if meta_lang is None:  # lang not specified
                return True
            elif meta_lang == self.lang.lang:  # matched to html_search_language
                return True

        return False

    def dispatch_visit(self, node: Node) -> None:
        if isinstance(node, nodes.comment):
            raise nodes.SkipNode
        elif isinstance(node, nodes.raw):
            if 'html' in node.get('format', '').split():
                # Some people might put content in raw HTML that should be searched,
                # so we just amateurishly strip HTML tags and index the remaining
                # content
                nodetext = re.sub(r'(?is)<style.*?</style>', '', node.astext())
                nodetext = re.sub(r'(?is)<script.*?</script>', '', nodetext)
                nodetext = re.sub(r'<[^<]+?>', '', nodetext)
                self.found_words.extend(self.lang.split(nodetext))
            raise nodes.SkipNode
        elif isinstance(node, nodes.Text):
            self.found_words.extend(self.lang.split(node.astext()))
        elif isinstance(node, nodes.title):
            self.found_title_words.extend(self.lang.split(node.astext()))
        elif isinstance(node, Element) and self.is_meta_keywords(node):
            keywords = node['content']
            keywords = [keyword.strip() for keyword in keywords.split(',')]
            self.found_words.extend(keywords)


class IndexBuilder:
    """
    Helper class that creates a search index based on the doctrees
    passed to the `feed` method.
    """
    formats = {
        'jsdump':   jsdump,
        'pickle':   pickle
    }

    def __init__(self, env: BuildEnvironment, lang: str, options: Dict, scoring: str) -> None:
        self.env = env
        self._titles: Dict[str, str] = {}           # docname -> title
        self._filenames: Dict[str, str] = {}        # docname -> filename
        self._mapping: Dict[str, Set[str]] = {}     # stemmed word -> set(docname)
        # stemmed words in titles -> set(docname)
        self._title_mapping: Dict[str, Set[str]] = {}
        self._stem_cache: Dict[str, str] = {}       # word -> stemmed word
        self._objtypes: Dict[Tuple[str, str], int] = {}     # objtype -> index
        # objtype index -> (domain, type, objname (localized))
        self._objnames: Dict[int, Tuple[str, str, str]] = {}
        # add language-specific SearchLanguage instance
        lang_class: Type[SearchLanguage] = languages.get(lang)

        # fallback; try again with language-code
        if lang_class is None and '_' in lang:
            lang_class = languages.get(lang.split('_')[0])

        if lang_class is None:
            self.lang: SearchLanguage = SearchEnglish(options)
        elif isinstance(lang_class, str):
            module, classname = lang_class.rsplit('.', 1)
            lang_class = getattr(import_module(module), classname)
            self.lang = lang_class(options)
        else:
            # it's directly a class (e.g. added by app.add_search_language)
            self.lang = lang_class(options)

        if scoring:
            with open(scoring, 'rb') as fp:
                self.js_scorer_code = fp.read().decode()
        else:
            self.js_scorer_code = ''
        self.js_splitter_code = splitter_code

    def load(self, stream: IO, format: Any) -> None:
        """Reconstruct from frozen data."""
        if isinstance(format, str):
            format = self.formats[format]
        frozen = format.load(stream)
        # if an old index is present, we treat it as not existing.
        if not isinstance(frozen, dict) or \
           frozen.get('envversion') != self.env.version:
            raise ValueError('old format')
        index2fn = frozen['docnames']
        self._filenames = dict(zip(index2fn, frozen['filenames']))
        self._titles = dict(zip(index2fn, frozen['titles']))

        def load_terms(mapping: Dict[str, Any]) -> Dict[str, Set[str]]:
            rv = {}
            for k, v in mapping.items():
                if isinstance(v, int):
                    rv[k] = {index2fn[v]}
                else:
                    rv[k] = {index2fn[i] for i in v}
            return rv

        self._mapping = load_terms(frozen['terms'])
        self._title_mapping = load_terms(frozen['titleterms'])
        # no need to load keywords/objtypes

    def dump(self, stream: IO, format: Any) -> None:
        """Dump the frozen index to a stream."""
        if isinstance(format, str):
            format = self.formats[format]
        format.dump(self.freeze(), stream)

    def get_objects(self, fn2index: Dict[str, int]
                    ) -> Dict[str, List[Tuple[int, int, int, str, str]]]:
        rv: Dict[str, List[Tuple[int, int, int, str, str]]] = {}
        otypes = self._objtypes
        onames = self._objnames
        for domainname, domain in sorted(self.env.domains.items()):
            for fullname, dispname, type, docname, anchor, prio in \
                    sorted(domain.get_objects()):
                if docname not in fn2index:
                    continue
                if prio < 0:
                    continue
                fullname = html.escape(fullname)
                dispname = html.escape(dispname)
                prefix, _, name = dispname.rpartition('.')
                plist = rv.setdefault(prefix, [])
                try:
                    typeindex = otypes[domainname, type]
                except KeyError:
                    typeindex = len(otypes)
                    otypes[domainname, type] = typeindex
                    otype = domain.object_types.get(type)
                    if otype:
                        # use str() to fire translation proxies
                        onames[typeindex] = (domainname, type,
                                             str(domain.get_type_name(otype)))
                    else:
                        onames[typeindex] = (domainname, type, type)
                if anchor == fullname:
                    shortanchor = ''
                elif anchor == type + '-' + fullname:
                    shortanchor = '-'
                else:
                    shortanchor = anchor
                plist.append((fn2index[docname], typeindex, prio, shortanchor, name))
        return rv

    def get_terms(self, fn2index: Dict) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
        rvs: Tuple[Dict[str, List[str]], Dict[str, List[str]]] = ({}, {})
        for rv, mapping in zip(rvs, (self._mapping, self._title_mapping)):
            for k, v in mapping.items():
                if len(v) == 1:
                    fn, = v
                    if fn in fn2index:
                        rv[k] = fn2index[fn]
                else:
                    rv[k] = sorted([fn2index[fn] for fn in v if fn in fn2index])
        return rvs

    def freeze(self) -> Dict[str, Any]:
        """Create a usable data structure for serializing."""
        docnames, titles = zip(*sorted(self._titles.items()))
        filenames = [self._filenames.get(docname) for docname in docnames]
        fn2index = {f: i for (i, f) in enumerate(docnames)}
        terms, title_terms = self.get_terms(fn2index)

        objects = self.get_objects(fn2index)  # populates _objtypes
        objtypes = {v: k[0] + ':' + k[1] for (k, v) in self._objtypes.items()}
        objnames = self._objnames
        return dict(docnames=docnames, filenames=filenames, titles=titles, terms=terms,
                    objects=objects, objtypes=objtypes, objnames=objnames,
                    titleterms=title_terms, envversion=self.env.version)

    def label(self) -> str:
        return "%s (code: %s)" % (self.lang.language_name, self.lang.lang)

    def prune(self, docnames: Iterable[str]) -> None:
        """Remove data for all docnames not in the list."""
        new_titles = {}
        new_filenames = {}
        for docname in docnames:
            if docname in self._titles:
                new_titles[docname] = self._titles[docname]
                new_filenames[docname] = self._filenames[docname]
        self._titles = new_titles
        self._filenames = new_filenames
        for wordnames in self._mapping.values():
            wordnames.intersection_update(docnames)
        for wordnames in self._title_mapping.values():
            wordnames.intersection_update(docnames)

    def feed(self, docname: str, filename: str, title: str, doctree: nodes.document) -> None:
        """Feed a doctree to the index."""
        self._titles[docname] = title
        self._filenames[docname] = filename

        visitor = WordCollector(doctree, self.lang)
        doctree.walk(visitor)

        # memoize self.lang.stem
        def stem(word: str) -> str:
            try:
                return self._stem_cache[word]
            except KeyError:
                self._stem_cache[word] = self.lang.stem(word).lower()
                return self._stem_cache[word]
        _filter = self.lang.word_filter

        for word in visitor.found_title_words:
            stemmed_word = stem(word)
            if _filter(stemmed_word):
                self._title_mapping.setdefault(stemmed_word, set()).add(docname)
            elif _filter(word): # stemmer must not remove words from search index
                self._title_mapping.setdefault(word, set()).add(docname)

        for word in visitor.found_words:
            stemmed_word = stem(word)
            # again, stemmer must not remove words from search index
            if not _filter(stemmed_word) and _filter(word):
                stemmed_word = word
            already_indexed = docname in self._title_mapping.get(stemmed_word, set())
            if _filter(stemmed_word) and not already_indexed:
                self._mapping.setdefault(stemmed_word, set()).add(docname)

    def context_for_searchtool(self) -> Dict[str, Any]:
        if self.lang.js_splitter_code:
            js_splitter_code = self.lang.js_splitter_code
        else:
            js_splitter_code = self.js_splitter_code

        return {
            'search_language_stemming_code': self.get_js_stemmer_code(),
            'search_language_stop_words': jsdump.dumps(sorted(self.lang.stopwords)),
            'search_scorer_tool': self.js_scorer_code,
            'search_word_splitter_code': js_splitter_code,
        }

    def get_js_stemmer_rawcodes(self) -> List[str]:
        """Returns a list of non-minified stemmer JS files to copy."""
        if self.lang.js_stemmer_rawcode:
            return [
                path.join(package_dir, 'search', 'non-minified-js', fname)
                for fname in ('base-stemmer.js', self.lang.js_stemmer_rawcode)
            ]
        else:
            return []

    def get_js_stemmer_rawcode(self) -> Optional[str]:
        return None

    def get_js_stemmer_code(self) -> str:
        """Returns JS code that will be inserted into language_data.js."""
        if self.lang.js_stemmer_rawcode:
            js_dir = path.join(package_dir, 'search', 'minified-js')
            with open(path.join(js_dir, 'base-stemmer.js')) as js_file:
                base_js = js_file.read()
            with open(path.join(js_dir, self.lang.js_stemmer_rawcode)) as js_file:
                language_js = js_file.read()
            return ('%s\n%s\nStemmer = %sStemmer;' %
                    (base_js, language_js, self.lang.language_name))
        else:
            return self.lang.js_stemmer_code
