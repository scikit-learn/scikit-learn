# Authors: Olivier Grisel <olivier.grisel@ensta.org>
#          Mathieu Blondel
#
# License: BSD Style.
"""Utilities to build feature vectors from text documents"""

from collections import defaultdict
import re
import unicodedata
import numpy as np
import scipy.sparse as sp
from ...base import BaseEstimator

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
    """Transform accentuated unicode symbols into their simple counterpart"""
    return ''.join((c for c in unicodedata.normalize('NFD', s)
                    if unicodedata.category(c) != 'Mn'))

def strip_tags(s):
    return re.compile(r"<([^>]+)>", flags=re.UNICODE).sub("", s)


class DefaultPreprocessor(object):

    def preprocess(self, text):
        return strip_accents(strip_tags(text.lower()))

    def __repr__(self):
        return "DefaultPreprocessor()"


DEFAULT_PREPROCESSOR = DefaultPreprocessor()


class WordNGramAnalyzer(BaseEstimator):
    """Simple analyzer: transform a text document into a sequence of word tokens

    This simple implementation does:
      - lower case conversion
      - unicode accents removal
      - token extraction using unicode regexp word bounderies for token of
        minimum size of 2 symbols
      - output token n-grams (unigram only by default)
    """

    token_pattern = re.compile(r"\b\w\w+\b", re.UNICODE)

    def __init__(self, charset='utf-8', min_n=1, max_n=1,
                 preprocessor=DEFAULT_PREPROCESSOR,
                 stop_words=ENGLISH_STOP_WORDS):
        self.charset = charset
        self.stop_words = stop_words
        self.min_n = min_n
        self.max_n = max_n
        self.preprocessor = preprocessor

    def analyze(self, text_document):
        if isinstance(text_document, file):
            text_document = text_document.read()

        if isinstance(text_document, str):
            text_document = text_document.decode(self.charset, 'ignore')

        text_document = self.preprocessor.preprocess(text_document)

        # word boundaries tokenizer
        tokens = self.token_pattern.findall(text_document)

        # handle token n-grams
        if self.min_n != 1 or self.max_n != 1:
            original_tokens = tokens
            tokens = []
            n_original_tokens = len(original_tokens)
            for n in xrange(self.min_n, self.max_n + 1):
                if n_original_tokens < n:
                    continue
                for i in xrange(n_original_tokens - n + 1):
                    tokens.append(" ".join(original_tokens[i: i + n]))

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


class BaseCountVectorizer(BaseEstimator):
    """Convert a collection of raw documents to a matrix of token counts

    This class can't be used directly, use either CountVectorizer or
    sparse.CountVectorizer.

    Parameters
    ----------
    analyzer: WordNGramAnalyzer or CharNGramAnalyzer, optional

    vocabulary: dict, optional
        A dictionary where keys are tokens and values are indices in the
        matrix.
        This is useful in order to fix the vocabulary in advance.

    dtype: type, optional
        Type of the matrix returned by fit_transform() or transform().
    """

    def __init__(self, analyzer=DEFAULT_ANALYZER, vocabulary={}, dtype=long):
        self.analyzer = analyzer
        self.vocabulary = vocabulary
        self.dtype = dtype

    def _init_matrix(self, shape):
        raise NotImplementedError

    def _build_vectors_and_vocab(self, raw_documents):
        vocab = {} # token => idx
        docs = []

        for doc in raw_documents:
            doc_dict = {} # idx => count

            for token in self.analyzer.analyze(doc):
                if not token in vocab:
                    vocab[token] = len(vocab)
                idx = vocab[token]
                doc_dict[idx] = doc_dict.get(idx, 0) + 1

            docs.append(doc_dict)

        # convert to a document-token matrix
        matrix = self._init_matrix((len(docs), len(vocab)))

        for i, doc_dict in enumerate(docs):
            for idx, count in doc_dict.iteritems():
                matrix[i, idx] = count

        return matrix, vocab

    def _build_vectors(self, raw_documents):
        # raw_documents is an iterable so we don't know its size in advance
        vectors = None

        for i, doc in enumerate(raw_documents):
            vector = self._init_matrix((1, len(self.vocabulary)))
            for token in self.analyzer.analyze(doc):
                try:
                    vector[0, self.vocabulary[token]] += 1
                except KeyError:
                    # ignore out-of-vocabulary tokens
                    pass

            if vectors is None:
                vectors = vector
            else:
                vstack = sp.vstack if sp.issparse(vector) else np.vstack
                vectors = vstack((vectors, vector))

        return vectors

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


class CountVectorizer(BaseCountVectorizer):
    """Convert a collection of raw documents to a matrix of token counts

    This implementation produces a dense representation of the counts using
    a numpy array.

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

    dtype: type, optional
        Type of the matrix returned by fit_transform() or transform().
    """


    def _init_matrix(self, shape):
        return np.zeros(shape, dtype=self.dtype)


class BaseTfidfTransformer(BaseEstimator):
    """Transform a count matrix to a TF or TF-IDF representation

    TF means term-frequency while TF-IDF means term-frequency times inverse
    document-frequency:

      http://en.wikipedia.org/wiki/TF-IDF

    This class can't be used directly, use either TfidfTransformer or
    sparse.TfidfTransformer.

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


class TfidfTransformer(BaseTfidfTransformer):
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

    def fit(self, X, y=None):
        """Learn the IDF vector (global term weights)

        Parameters
        ----------
        X: array, [n_samples, n_features]
            a matrix of term/token counts

        """

        if self.use_idf:
            # how many documents include each token?
            idc = np.sum(X > 0, axis=0)
            self.idf = np.log(float(X.shape[0]) / idc)

        return self

    def transform(self, X, copy=True):
        """Transform a count matrix to a TF or TF-IDF representation

        Parameters
        ----------
        X: array, [n_samples, n_features]
            a matrix of term/token counts

        Returns
        -------
        vectors: array, [n_samples, n_features]
        """
        X = np.array(X, dtype=np.float64, copy=copy)

        if self.use_tf:
            # term-frequencies (normalized counts)
            X /= np.sum(X, axis=1)[:,np.newaxis]

        if self.use_idf:
            X *= self.idf

        return X


class BaseVectorizer(BaseEstimator):
    """Convert a collection of raw documents to a matrix

    This class can't be used directly, use either Vectorizer or
    sparse.Vectorizer.
    """

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


class Vectorizer(BaseVectorizer):
    """Convert a collection of raw documents to a matrix

    Equivalent to CountVectorizer followed by TfidfTransformer.
    """

    def __init__(self,
                 analyzer=DEFAULT_ANALYZER,
                 use_tf=True,
                 use_idf=True):
        self.tc = CountVectorizer(analyzer, dtype=np.float64)
        self.tfidf = TfidfTransformer(use_tf, use_idf)

