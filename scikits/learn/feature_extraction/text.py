# Authors: Olivier Grisel <olivier.grisel@ensta.org>
#          Mathieu Blondel
#
# License: BSD Style.
"""Utilities to build dense feature vectors from text documents"""

from operator import itemgetter
import re
import unicodedata
import numpy as np
from ..base import BaseEstimator
from ..preprocessing.sparse import Normalizer

ENGLISH_STOP_WORDS = set([
    "a", "about", "above", "across", "after", "afterwards", "again", "against",
    "all", "almost", "alone", "along", "already", "also", "although", "always",
    "am", "among", "amongst", "amoungst", "amount", "an", "and", "another",
    "any", "anyhow", "anyone", "anything", "anyway", "anywhere", "are",
    "around", "as", "at", "back", "be", "became", "because", "become",
    "becomes", "becoming", "been", "before", "beforehand", "behind", "being",
    "below", "beside", "besides", "between", "beyond", "bill", "both", "bottom",
    "but", "by", "call", "can", "cannot", "cant", "co", "computer", "con",
    "could", "couldnt", "cry", "de", "describe", "detail", "do", "done", "down",
    "due", "during", "each", "eg", "eight", "either", "eleven", "else",
    "elsewhere", "empty", "enough", "etc", "even", "ever", "every", "everyone",
    "everything", "everywhere", "except", "few", "fifteen", "fify", "fill",
    "find", "fire", "first", "five", "for", "former", "formerly", "forty",
    "found", "four", "from", "front", "full", "further", "get", "give", "go",
    "had", "has", "hasnt", "have", "he", "hence", "her", "here", "hereafter",
    "hereby", "herein", "hereupon", "hers", "herself", "him", "himself", "his",
    "how", "however", "hundred", "i", "ie", "if", "in", "inc", "indeed",
    "interest", "into", "is", "it", "its", "itself", "keep", "last", "latter",
    "latterly", "least", "less", "ltd", "made", "many", "may", "me",
    "meanwhile", "might", "mill", "mine", "more", "moreover", "most", "mostly",
    "move", "much", "must", "my", "myself", "name", "namely", "neither", "never",
    "nevertheless", "next", "nine", "no", "nobody", "none", "noone", "nor",
    "not", "nothing", "now", "nowhere", "of", "off", "often", "on", "once",
    "one", "only", "onto", "or", "other", "others", "otherwise", "our", "ours",
    "ourselves", "out", "over", "own", "part", "per", "perhaps", "please",
    "put", "rather", "re", "same", "see", "seem", "seemed", "seeming", "seems",
    "serious", "several", "she", "should", "show", "side", "since", "sincere",
    "six", "sixty", "so", "some", "somehow", "someone", "something", "sometime",
    "sometimes", "somewhere", "still", "such", "system", "take", "ten", "than",
    "that", "the", "their", "them", "themselves", "then", "thence", "there",
    "thereafter", "thereby", "therefore", "therein", "thereupon", "these",
    "they", "thick", "thin", "third", "this", "those", "though", "three",
    "through", "throughout", "thru", "thus", "to", "together", "too", "top",
    "toward", "towards", "twelve", "twenty", "two", "un", "under", "until",
    "up", "upon", "us", "very", "via", "was", "we", "well", "were", "what",
    "whatever", "when", "whence", "whenever", "where", "whereafter", "whereas",
    "whereby", "wherein", "whereupon", "wherever", "whether", "which", "while",
    "whither", "who", "whoever", "whole", "whom", "whose", "why", "will",
    "with", "within", "without", "would", "yet", "you", "your", "yours",
    "yourself", "yourselves"])


def strip_accents(s):
    """Transform accentuated unicode symbols into their simple counterpart

    Warning: the python-level loop and join operations make this implementation
    20 times slower than the to_ascii basic normalization.
    """
    return u''.join([c for c in unicodedata.normalize('NFKD', s)
                     if not unicodedata.combining(c)])


def to_ascii(s):
    """Transform accentuated unicode symbols into ascii or nothing

    Warning: this solution is only suited for roman languages that have a direct
    transliteration to ASCII symbols.

    A better solution would be to use transliteration based on a precomputed
    unidecode map to be used by translate as explained here:

        http://stackoverflow.com/questions/2854230/

    """
    nkfd_form = unicodedata.normalize('NFKD', s)
    only_ascii = nkfd_form.encode('ASCII', 'ignore')
    return only_ascii


def strip_tags(s):
    return re.compile(r"<([^>]+)>", flags=re.UNICODE).sub("", s)


class RomanPreprocessor(object):
    """Fast preprocessor suitable for roman languages"""

    def preprocess(self, unicode_text):
        return to_ascii(strip_tags(unicode_text.lower()))

    def __repr__(self):
        return "RomanPreprocessor()"


DEFAULT_PREPROCESSOR = RomanPreprocessor()

DEFAULT_TOKEN_PATTERN = r"\b\w\w+\b"

class WordNGramAnalyzer(BaseEstimator):
    """Simple analyzer: transform a text document into a sequence of word tokens

    This simple implementation does:
      - lower case conversion
      - unicode accents removal
      - token extraction using unicode regexp word bounderies for token of
        minimum size of 2 symbols (by default)
      - output token n-grams (unigram only by default)
    """

    def __init__(self, charset='utf-8', min_n=1, max_n=1,
                 preprocessor=DEFAULT_PREPROCESSOR,
                 stop_words=ENGLISH_STOP_WORDS,
                 token_pattern=DEFAULT_TOKEN_PATTERN):
        self.charset = charset
        self.stop_words = stop_words
        self.min_n = min_n
        self.max_n = max_n
        self.preprocessor = preprocessor
        self.token_pattern = token_pattern

    def analyze(self, text_document):
        if hasattr(text_document, 'read'):
            # ducktype for file-like objects
            text_document = text_document.read()

        if isinstance(text_document, str):
            text_document = text_document.decode(self.charset, 'ignore')

        text_document = self.preprocessor.preprocess(text_document)

        # word boundaries tokenizer (cannot compile it in the __init__ because
        # we want support for pickling and runtime parameter fitting)
        compiled = re.compile(self.token_pattern, re.UNICODE)
        tokens = compiled.findall(text_document)

        # handle token n-grams
        if self.min_n != 1 or self.max_n != 1:
            original_tokens = tokens
            tokens = []
            n_original_tokens = len(original_tokens)
            for n in xrange(self.min_n, self.max_n + 1):
                if n_original_tokens < n:
                    continue
                for i in xrange(n_original_tokens - n + 1):
                    tokens.append(u" ".join(original_tokens[i: i + n]))

        # handle stop words
        if self.stop_words is not None:
            tokens = [w for w in tokens if w not in self.stop_words]

        return tokens


class CharNGramAnalyzer(BaseEstimator):
    """Compute character n-grams features of a text document

    This analyzer is interesting since it is language agnostic and will work
    well even for language where word segmentation is not as trivial as English
    such as Chinese and German for instance.

    Because of this, it can be considered a basic morphological analyzer.
    """

    white_spaces = re.compile(r"\s\s+")

    def __init__(self, charset='utf-8', preprocessor=DEFAULT_PREPROCESSOR,
                 min_n=3, max_n=6):
        self.charset = charset
        self.min_n = min_n
        self.max_n = max_n
        self.preprocessor = preprocessor

    def analyze(self, text_document):
        if hasattr(text_document, 'read'):
            # ducktype for file-like objects
            text_document = text_document.read()

        if isinstance(text_document, str):
            text_document = text_document.decode(self.charset, 'ignore')

        text_document = self.preprocessor.preprocess(text_document)

        # normalize white spaces
        text_document = self.white_spaces.sub(" ", text_document)

        text_len = len(text_document)
        ngrams = []
        for n in xrange(self.min_n, self.max_n + 1):
            if text_len < n:
                continue
            for i in xrange(text_len - n):
                ngrams.append(text_document[i: i + n])
        return ngrams


DEFAULT_ANALYZER = WordNGramAnalyzer(min_n=1, max_n=1)


class CountVectorizer(BaseEstimator):
    """Convert a collection of raw documents to a matrix of token counts

    This implementation produces a sparse representation of the counts using
    scipy.sparse.coo_matrix.

    If you do not provide an a-priori dictionary and you do not use
    an analyzer that does some kind of feature selection then the number of
    features (the vocabulary size found by analysing the data) might be very
    large and the count vectors might not fit in memory.

    For this case it is either recommended to use the sparse.CountVectorizer
    variant of this class or a HashingVectorizer that will reduce the
    dimensionality to an arbitrary number by using random projection.

    Parameters
    ----------
    analyzer: WordNGramAnalyzer or CharNGramAnalyzer, optional

    vocabulary: dict, optional
        A dictionary where keys are tokens and values are indices in the
        matrix.

        This is useful in order to fix the vocabulary in advance.

    max_df : float in range [0.0, 1.0], optional, 0.5 by default
        When building the vocabulary ignore terms that have a term frequency
        high than the given threshold (corpus specific stop words).

        This parameter is ignored if vocabulary is not None.

    max_features : optional, None by default
        If not None, build a vocabulary that only consider the top
        max_features ordered by term frequency across the corpus.

        This parameter is ignored if vocabulary is not None.

    dtype: type, optional
        Type of the matrix returned by fit_transform() or transform().
    """

    def __init__(self, analyzer=DEFAULT_ANALYZER, vocabulary={}, max_df=1.0,
                 max_features=None, dtype=long):
        self.analyzer = analyzer
        self.vocabulary = vocabulary
        self.dtype = dtype
        self.max_df = max_df
        self.max_features = max_features

    def _term_count_dicts_to_matrix(self, term_count_dicts, vocabulary):

        import scipy.sparse as sp
        i_indices = []
        j_indices = []
        values = []

        for i, term_count_dict in enumerate(term_count_dicts):
            for term, count in term_count_dict.iteritems():
                j = vocabulary.get(term)
                if j is not None:
                    i_indices.append(i)
                    j_indices.append(j)
                    values.append(count)
            # free memory as we go
            term_count_dict.clear()

        shape = (len(term_count_dicts), max(vocabulary.itervalues()) + 1)
        return sp.coo_matrix((values, (i_indices, j_indices)),
                             shape=shape, dtype=self.dtype)


    def _build_vectors_and_vocab(self, raw_documents):
        """Analyze documents, build vocabulary and vectorize"""

        # result of document conversion to term_count_dict
        term_counts_per_doc = []
        term_counts = {}

        # term counts across entire corpus (count each term maximum once per
        # document)
        document_counts = {}

        max_df = self.max_df
        max_features = self.max_features

        # TODO: parallelize the following loop with joblib
        for doc in raw_documents:
            term_count_dict = {} # term => count in doc

            for term in self.analyzer.analyze(doc):
                term_count_dict[term] = term_count_dict.get(term, 0) + 1
                term_counts[term] = term_counts.get(term, 0) + 1

            if max_df is not None:
                for term in term_count_dict.iterkeys():
                    document_counts[term] = document_counts.get(term, 0) + 1

            term_counts_per_doc.append(term_count_dict)

        n_doc = len(term_counts_per_doc)

        # filter out stop words: terms that occur in almost all documents
        stop_words = set()
        if max_df is not None:
            max_document_count = max_df * n_doc
            for t, dc in sorted(document_counts.iteritems(), key=itemgetter(1),
                                reverse=True):
                if dc < max_document_count:
                    break
                stop_words.add(t)

        # list the terms that should be part of the vocabulary
        if max_features is not None:
            # extract the most frequent terms for the vocabulary
            terms = set()
            for t, tc in sorted(term_counts.iteritems(), key=itemgetter(1),
                                reverse=True):
                if t not in stop_words:
                    terms.add(t)
                if len(terms) >= max_features:
                    break
        else:
            terms = set(term_counts.keys())
            terms -= stop_words

        # convert to a document-token matrix
        vocabulary = dict(((t, i) for i, t in enumerate(terms))) # token => idx

        # the term_counts and document_counts might be useful statistics, are
        # we really sure want we want to drop them? They take some memory but
        # can be useful for corpus introspection

        matrix = self._term_count_dicts_to_matrix(term_counts_per_doc, vocabulary)
        return matrix, vocabulary

    def _build_vectors(self, raw_documents):
        """Analyze documents and vectorize using existing vocabulary"""
        # raw_documents is an iterable so we don't know its size in advance

        # result of document conversion to term_count_dict
        term_counts_per_doc = []

        # TODO: parallelize the following loop with joblib
        for doc in raw_documents:
            term_count_dict = {} # term => count in doc

            for term in self.analyzer.analyze(doc):
                term_count_dict[term] = term_count_dict.get(term, 0) + 1

            term_counts_per_doc.append(term_count_dict)

        # now that we know the document we can allocate the vectors matrix at
        # once and fill it with the term counts collected as a temporary list
        # of dict
        return self._term_count_dicts_to_matrix(
            term_counts_per_doc, self.vocabulary)

    def fit(self, raw_documents, y=None):
        """Learn a vocabulary dictionary of all tokens in the raw documents

        Parameters
        ----------

        raw_documents: iterable
            an iterable which yields either str, unicode or file objects

        Returns
        -------
        self
        """
        self.fit_transform(raw_documents)
        return self

    def fit_transform(self, raw_documents, y=None):
        """Learn the vocabulary dictionary and return the count vectors

        This is more efficient than calling fit followed by transform.

        Parameters
        ----------

        raw_documents: iterable
            an iterable which yields either str, unicode or file objects

        Returns
        -------
        vectors: array, [n_samples, n_features]
        """
        vectors, self.vocabulary = self._build_vectors_and_vocab(raw_documents)
        return vectors

    def transform(self, raw_documents):
        """Extract token counts out of raw text documents

        Parameters
        ----------

        raw_documents: iterable
            an iterable which yields either str, unicode or file objects

        Returns
        -------
        vectors: array, [n_samples, n_features]
        """
        if len(self.vocabulary) == 0:
            raise ValueError, "No vocabulary dictionary available..."

        return self._build_vectors(raw_documents)


class TfidfTransformer(BaseEstimator):
    """Transform a count matrix to a TF or TF-IDF representation

    TF means term-frequency while TF-IDF means term-frequency times inverse
    document-frequency:

      http://en.wikipedia.org/wiki/TF-IDF

    The goal of using TF-IDF instead of the raw frequencies of occurrence of a
    token in a given document is to scale down the impact of tokens that occur
    very frequently in a given corpus and that are hence empirically less
    informative than feature that occur in a small fraction of the training
    corpus.

    TF-IDF can be seen as a smooth alternative to the stop words filtering.

    Parameters
    ----------

    use_tf: boolean
        enable term-frequency normalization

    use_idf: boolean
        enable inverse-document-frequency reweighting
    """

    def __init__(self, use_tf=True, use_idf=True):
        self.use_tf = use_tf
        self.use_idf = use_idf
        self.idf = None

    def fit(self, X, y=None):
        """Learn the IDF vector (global term weights)

        Parameters
        ----------
        X: sparse matrix, [n_samples, n_features]
            a matrix of term/token counts

        """
        n_samples, n_features = X.shape
        if self.use_idf:
            # how many documents include each token?
            idc = np.zeros(n_features, dtype=np.float64)
            for doc, token in zip(*X.nonzero()):
                idc[token] += 1
            self.idf = np.log(float(X.shape[0]) / idc)

        return self

    def transform(self, X, copy=True):
        """Transform a count matrix to a TF or TF-IDF representation

        Parameters
        ----------
        X: sparse matrix, [n_samples, n_features]
            a matrix of term/token counts

        Returns
        -------
        vectors: sparse matrix, [n_samples, n_features]
        """
        import scipy.sparse as sp
        X = sp.csr_matrix(X, dtype=np.float64, copy=copy)
        n_samples, n_features = X.shape

        if self.use_tf:
            X = Normalizer().transform(X)

        if self.use_idf:
            d = sp.lil_matrix((len(self.idf), len(self.idf)))
            d.setdiag(self.idf)
            # *= doesn't work
            X = X * d

        return X


class Vectorizer(BaseEstimator):
    """Convert a collection of raw documents to a matrix

    Equivalent to CountVectorizer followed by TfidfTransformer.
    """

    def __init__(self, analyzer=DEFAULT_ANALYZER, max_df=1.0,
                 max_features=None, use_tf=True, use_idf=True):
        self.tc = CountVectorizer(analyzer, max_df=max_df,
                                  max_features=max_features,
                                  dtype=np.float64)
        self.tfidf = TfidfTransformer(use_tf, use_idf)

    def fit(self, raw_documents):
        X = self.tc.fit_transform(raw_documents)
        self.tfidf.fit(X)
        return self

    def fit_transform(self, raw_documents):
        """
        Learn the representation and return the vectors.

        Parameters
        ----------

        raw_documents: iterable
            an iterable which yields either str, unicode or file objects

        Returns
        -------
        vectors: array, [n_samples, n_features]
        """
        X = self.tc.fit_transform(raw_documents)
        # X is already a transformed view of raw_documents so
        # we set copy to False
        return self.tfidf.fit(X).transform(X, copy=False)

    def transform(self, raw_documents, copy=True):
        """
        Return the vectors.

        Parameters
        ----------

        raw_documents: iterable
            an iterable which yields either str, unicode or file objects

        Returns
        -------
        vectors: array, [n_samples, n_features]
        """
        X = self.tc.transform(raw_documents)
        return self.tfidf.transform(X, copy)

    def _get_vocab(self):
        return self.tc.vocabulary

    vocabulary = property(_get_vocab)
