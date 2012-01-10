# -*- coding: utf-8 -*-
# Authors: Olivier Grisel <olivier.grisel@ensta.org>
#          Mathieu Blondel <mathieu@mblondel.org>
#          Lars Buitinck <L.J.Buitinck@uva.nl>
#
# License: BSD Style.
"""
The :mod:`sklearn.feature_extraction.text` submodule gathers utilities to
build feature vectors from text documents.
"""

import re
import unicodedata
import numpy as np
import scipy.sparse as sp
from ..base import BaseEstimator, TransformerMixin
from ..preprocessing import normalize
from ..utils.fixes import Counter

# This list of English stop words is taken from the "Glasgow Information
# Retrieval Group". The original list can be found at
# http://ir.dcs.gla.ac.uk/resources/linguistic_utils/stop_words
ENGLISH_STOP_WORDS = frozenset([
    "a", "about", "above", "across", "after", "afterwards", "again", "against",
    "all", "almost", "alone", "along", "already", "also", "although", "always",
    "am", "among", "amongst", "amoungst", "amount", "an", "and", "another",
    "any", "anyhow", "anyone", "anything", "anyway", "anywhere", "are",
    "around", "as", "at", "back", "be", "became", "because", "become",
    "becomes", "becoming", "been", "before", "beforehand", "behind", "being",
    "below", "beside", "besides", "between", "beyond", "bill", "both",
    "bottom", "but", "by", "call", "can", "cannot", "cant", "co", "con",
    "could", "couldnt", "cry", "de", "describe", "detail", "do", "done",
    "down", "due", "during", "each", "eg", "eight", "either", "eleven", "else",
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
    "move", "much", "must", "my", "myself", "name", "namely", "neither",
    "never", "nevertheless", "next", "nine", "no", "nobody", "none", "noone",
    "nor", "not", "nothing", "now", "nowhere", "of", "off", "often", "on",
    "once", "one", "only", "onto", "or", "other", "others", "otherwise", "our",
    "ours", "ourselves", "out", "over", "own", "part", "per", "perhaps",
    "please", "put", "rather", "re", "same", "see", "seem", "seemed",
    "seeming", "seems", "serious", "several", "she", "should", "show", "side",
    "since", "sincere", "six", "sixty", "so", "some", "somehow", "someone",
    "something", "sometime", "sometimes", "somewhere", "still", "such",
    "system", "take", "ten", "than", "that", "the", "their", "them",
    "themselves", "then", "thence", "there", "thereafter", "thereby",
    "therefore", "therein", "thereupon", "these", "they", "thick", "thin",
    "third", "this", "those", "though", "three", "through", "throughout",
    "thru", "thus", "to", "together", "too", "top", "toward", "towards",
    "twelve", "twenty", "two", "un", "under", "until", "up", "upon", "us",
    "very", "via", "was", "we", "well", "were", "what", "whatever", "when",
    "whence", "whenever", "where", "whereafter", "whereas", "whereby",
    "wherein", "whereupon", "wherever", "whether", "which", "while", "whither",
    "who", "whoever", "whole", "whom", "whose", "why", "will", "with",
    "within", "without", "would", "yet", "you", "your", "yours", "yourself",
    "yourselves"])


def strip_accents(s):
    """Transform accentuated unicode symbols into their simple counterpart

    Warning: the python-level loop and join operations make this implementation
    20 times slower than the to_ascii basic normalization.
    """
    return u''.join([c for c in unicodedata.normalize('NFKD', s)
                     if not unicodedata.combining(c)])


def to_ascii(s):
    """Transform accentuated unicode symbols into ascii or nothing

    Warning: this solution is only suited for languages that have a direct
    transliteration to ASCII symbols.

    A better solution would be to use transliteration based on a precomputed
    unidecode map to be used by translate as explained here:

        http://stackoverflow.com/questions/2854230/

    """
    nkfd_form = unicodedata.normalize('NFKD', s)
    only_ascii = nkfd_form.encode('ASCII', 'ignore')
    return only_ascii


def strip_tags(s):
    return re.compile(ur"<([^>]+)>", flags=re.UNICODE).sub(u"", s)


class RomanPreprocessor(object):
    """Fast preprocessor suitable for Latin alphabet text"""

    def preprocess(self, unicode_text):
        """Preprocess strings"""
        return to_ascii(strip_tags(unicode_text.lower()))

    def __repr__(self):
        return "RomanPreprocessor()"


DEFAULT_PREPROCESSOR = RomanPreprocessor()

DEFAULT_TOKEN_PATTERN = ur"\b\w\w+\b"


def _check_stop_list(stop):
    if stop == "english":
        return ENGLISH_STOP_WORDS
    elif isinstance(stop, str) or isinstance(stop, unicode):
        raise ValueError("not a built-in stop list: %s" % stop)
    else:               # assume it's a collection
        return stop


class WordNGramAnalyzer(BaseEstimator):
    """Simple analyzer: transform text document into a sequence of word tokens

    This simple implementation does:
      - lower case conversion
      - unicode accents removal
      - token extraction using unicode regexp word bounderies for token of
        minimum size of 2 symbols (by default)
      - output token n-grams (unigram only by default)

    The stop words argument may be "english" for a built-in list of English
    stop words or a collection of strings. Note that stop word filtering is
    performed after preprocessing, which may include accent stripping.
    """

    def __init__(self, charset='utf-8', min_n=1, max_n=1,
                 preprocessor=DEFAULT_PREPROCESSOR,
                 stop_words="english",
                 token_pattern=DEFAULT_TOKEN_PATTERN):
        self.charset = charset
        self.stop_words = stop_words
        self.min_n = min_n
        self.max_n = max_n
        self.preprocessor = preprocessor
        self.token_pattern = token_pattern

    def analyze(self, text_document):
        """From documents to token"""
        if hasattr(text_document, 'read'):
            # ducktype for file-like objects
            text_document = text_document.read()

        if isinstance(text_document, bytes):
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
            for n in xrange(self.min_n,
                            min(self.max_n + 1, n_original_tokens + 1)):
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

    white_spaces = re.compile(ur"\s\s+")

    def __init__(self, charset='utf-8', preprocessor=DEFAULT_PREPROCESSOR,
                 min_n=3, max_n=6):
        self.charset = charset
        self.min_n = min_n
        self.max_n = max_n
        self.preprocessor = preprocessor

    def analyze(self, text_document):
        """From documents to token"""
        if hasattr(text_document, 'read'):
            # ducktype for file-like objects
            text_document = text_document.read()

        if isinstance(text_document, bytes):
            text_document = text_document.decode(self.charset, 'ignore')

        text_document = self.preprocessor.preprocess(text_document)

        # normalize white spaces
        text_document = self.white_spaces.sub(u" ", text_document)

        text_len = len(text_document)
        ngrams = []
        for n in xrange(self.min_n, min(self.max_n + 1, text_len + 1)):
            for i in xrange(text_len - n + 1):
                ngrams.append(text_document[i: i + n])
        return ngrams


DEFAULT_ANALYZER = WordNGramAnalyzer(min_n=1, max_n=1)


class CountVectorizer(BaseEstimator):
    """Convert a collection of raw documents to a matrix of token counts

    This implementation produces a sparse representation of the counts using
    scipy.sparse.coo_matrix.

    If you do not provide an a-priori dictionary and you do not use an analyzer
    that does some kind of feature selection then the number of features will
    be equal to the vocabulary size found by analysing the data. The default
    analyzer does simple stop word filtering for English.

    Parameters
    ----------
    analyzer: WordNGramAnalyzer or CharNGramAnalyzer, optional

    vocabulary: dict or iterable, optional
        Either a dictionary where keys are tokens and values are indices in
        the matrix, or an iterable over terms (in which case the indices are
        determined by the iteration order as per enumerate).

        This is useful in order to fix the vocabulary in advance.

    max_df : float in range [0.0, 1.0], optional, 1.0 by default
        When building the vocabulary ignore terms that have a term frequency
        strictly higher than the given threshold (corpus specific stop words).

        This parameter is ignored if vocabulary is not None.

    max_features : optional, None by default
        If not None, build a vocabulary that only consider the top
        max_features ordered by term frequency across the corpus.

        This parameter is ignored if vocabulary is not None.

    dtype: type, optional
        Type of the matrix returned by fit_transform() or transform().
    """

    def __init__(self, analyzer=DEFAULT_ANALYZER, vocabulary=None, max_df=1.0,
                 max_features=None, dtype=long):
        self.analyzer = analyzer
        self.fit_vocabulary = vocabulary is None
        if vocabulary is not None and not isinstance(vocabulary, dict):
            vocabulary = dict((t, i) for i, t in enumerate(vocabulary))
        self.vocabulary = vocabulary
        self.dtype = dtype
        self.max_df = max_df
        self.max_features = max_features

    def _term_count_dicts_to_matrix(self, term_count_dicts):
        i_indices = []
        j_indices = []
        values = []
        vocabulary = self.vocabulary

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
        if not self.fit_vocabulary:
            return self.transform(raw_documents)

        # result of document conversion to term count dicts
        term_counts_per_doc = []
        term_counts = Counter()

        # term counts across entire corpus (count each term maximum once per
        # document)
        document_counts = Counter()

        max_df = self.max_df
        max_features = self.max_features

        # TODO: parallelize the following loop with joblib?
        # (see XXX up ahead)
        for doc in raw_documents:
            term_count_current = Counter(self.analyzer.analyze(doc))
            term_counts.update(term_count_current)

            if max_df < 1.0:
                document_counts.update(term_count_current.iterkeys())

            term_counts_per_doc.append(term_count_current)

        n_doc = len(term_counts_per_doc)

        # filter out stop words: terms that occur in almost all documents
        if max_df < 1.0:
            max_document_count = max_df * n_doc
            stop_words = set(t for t, dc in document_counts.iteritems()
                               if dc > max_document_count)
        else:
            stop_words = set()

        # list the terms that should be part of the vocabulary
        if max_features is None:
            terms = set(term_counts) - stop_words
        else:
            # extract the most frequent terms for the vocabulary
            terms = set()
            for t, tc in term_counts.most_common():
                if t not in stop_words:
                    terms.add(t)
                if len(terms) >= max_features:
                    break

        # convert to a document-token matrix
        self.vocabulary = dict(((t, i) for i, t in enumerate(terms)))

        # the term_counts and document_counts might be useful statistics, are
        # we really sure want we want to drop them? They take some memory but
        # can be useful for corpus introspection

        return self._term_count_dicts_to_matrix(term_counts_per_doc)

    def transform(self, raw_documents):
        """Extract token counts out of raw text documents using the vocabulary
        fitted with fit or the one provided in the constructor.

        Parameters
        ----------
        raw_documents: iterable
            an iterable which yields either str, unicode or file objects

        Returns
        -------
        vectors: sparse matrix, [n_samples, n_features]
        """
        if not self.vocabulary:
            raise ValueError("Vocabulary wasn't fitted or is empty!")

        # raw_documents is an iterable so we don't know its size in advance

        # XXX @larsmans tried to parallelize the following loop with joblib.
        # The result was some 20% slower than the serial version.
        term_counts_per_doc = [Counter(self.analyzer.analyze(doc))
                               for doc in raw_documents]
        return self._term_count_dicts_to_matrix(term_counts_per_doc)

    def inverse_transform(self, X):
        """Return terms per document with nonzero entries in X.

        Parameters
        ----------
        X : {array, sparse matrix}, shape = [n_samples, n_features]

        Returns
        -------
        X_inv : list of arrays, len = n_samples
            List of arrays of terms.
        """
        if sp.isspmatrix_coo(X):  # COO matrix is not indexable
            X = X.tocsr()
        elif not sp.issparse(X):
            # We need to convert X to a matrix, so that the indexing
            # returns 2D objects
            X = np.asmatrix(X)
        n_samples = X.shape[0]

        terms = np.array(self.vocabulary.keys())
        indices = np.array(self.vocabulary.values())
        inverse_vocabulary = terms[np.argsort(indices)]

        return [inverse_vocabulary[X[i, :].nonzero()[1]].ravel()
                for i in xrange(n_samples)]


class TfidfTransformer(BaseEstimator, TransformerMixin):
    """Transform a count matrix to a normalized tf or tf–idf representation

    Tf means term-frequency while tf–idf means term-frequency times inverse
    document-frequency. This is a common term weighting scheme in information
    retrieval, that has also found good use in document classification.

    The goal of using tf–idf instead of the raw frequencies of occurrence of a
    token in a given document is to scale down the impact of tokens that occur
    very frequently in a given corpus and that are hence empirically less
    informative than features that occur in a small fraction of the training
    corpus.

    In the SMART notation used in IR, this class implements several tf–idf
    variants. Tf is always "n" (natural), idf is "t" iff use_idf is given,
    "n" otherwise, and normalization is "c" iff norm='l2', "n" iff norm=None.

    Parameters
    ----------
    norm : 'l1', 'l2' or None, optional
        Norm used to normalize term vectors. None for no normalization.

    use_idf : boolean, optional
        Enable inverse-document-frequency reweighting.

    smooth_idf : boolean, optional
        Smooth idf weights by adding one to document frequencies, as if an
        extra document was seen containing every term in the collection
        exactly once. Prevents zero divisions.

    Notes
    -----
    **References**:

    .. [Yates2011] `R. Baeza-Yates and B. Ribeiro-Neto (2011). Modern
                   Information Retrieval. Addison Wesley, pp. 68–74.`

    .. [MSR2008] `C.D. Manning, H. Schütze and P. Raghavan (2008). Introduction
                 to Information Retrieval. Cambridge University Press,
                 pp. 121–125.`
    """

    def __init__(self, norm='l2', use_idf=True, smooth_idf=True):
        self.norm = norm
        self.use_idf = use_idf
        self.smooth_idf = smooth_idf
        self.idf_ = None

    def fit(self, X, y=None):
        """Learn the idf vector (global term weights)

        Parameters
        ----------
        X: sparse matrix, [n_samples, n_features]
            a matrix of term/token counts
        """
        if self.use_idf:
            n_samples, n_features = X.shape
            df = np.bincount(X.nonzero()[1])
            if df.shape[0] < n_features:
                # bincount might return fewer bins than there are features
                df = np.concatenate([df, np.zeros(n_features - df.shape[0])])
            df += int(self.smooth_idf)
            self.idf_ = np.log(float(n_samples) / df)

        return self

    def transform(self, X, copy=True):
        """Transform a count matrix to a tf or tf–idf representation

        Parameters
        ----------
        X: sparse matrix, [n_samples, n_features]
            a matrix of term/token counts

        Returns
        -------
        vectors: sparse matrix, [n_samples, n_features]
        """
        X = sp.csr_matrix(X, dtype=np.float64, copy=copy)
        n_samples, n_features = X.shape

        if self.use_idf:
            expected_n_features = self.idf_.shape[0]
            if n_features != expected_n_features:
                raise ValueError("Input has n_features=%d while the model"
                                 " has been trained with n_features=%d" % (
                                     n_features, expected_n_features))
            d = sp.lil_matrix((n_features, n_features))
            d.setdiag(self.idf_)
            # *= doesn't work
            X = X * d

        if self.norm:
            X = normalize(X, norm=self.norm, copy=False)

        return X


class Vectorizer(BaseEstimator):
    """Convert a collection of raw documents to a matrix

    Equivalent to CountVectorizer followed by TfidfTransformer.
    """

    def __init__(self, analyzer=DEFAULT_ANALYZER, max_df=1.0,
                 max_features=None, norm='l2', use_idf=True, smooth_idf=True):
        self.tc = CountVectorizer(analyzer, max_df=max_df,
                                  max_features=max_features,
                                  dtype=np.float64)
        self.tfidf = TfidfTransformer(norm=norm, use_idf=use_idf,
                                      smooth_idf=smooth_idf)

    def fit(self, raw_documents):
        """Learn a conversion law from documents to array data"""
        X = self.tc.fit_transform(raw_documents)
        self.tfidf.fit(X)
        return self

    def fit_transform(self, raw_documents, y=None):
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
        """Transform raw text documents to tf–idf vectors

        Parameters
        ----------
        raw_documents: iterable
            an iterable which yields either str, unicode or file objects

        Returns
        -------
        vectors: sparse matrix, [n_samples, n_features]
        """
        X = self.tc.transform(raw_documents)
        return self.tfidf.transform(X, copy)

    def inverse_transform(self, X):
        """Return terms per document with nonzero entries in X.

        Parameters
        ----------
        X : {array, sparse matrix}, shape = [n_samples, n_features]

        Returns
        -------
        X_inv : list of arrays, len = n_samples
            List of arrays of terms.
        """
        return self.tc.inverse_transform(X)

    vocabulary = property(lambda self: self.tc.vocabulary)
    analyzer = property(lambda self: self.tc.analyzer)
