# Author: Olivier Grisel <olivier.grisel@ensta.org>
#
# License: BSD Style.
"""Utilities to build feature vectors from text documents"""

from collections import defaultdict
import re
import unicodedata
import numpy as np
import scipy.sparse as sp
from ..base import BaseEstimator
from ..pipeline import Pipeline


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


DEFAULT_FILTERS = [lambda x: x.lower(), strip_tags, strip_accents]


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
                 filters=DEFAULT_FILTERS, stop_words=None):
        self.charset = charset
        self.stop_words = stop_words
        self.min_n = min_n
        self.max_n = max_n
        self.filters = filters

    def analyze(self, text_document):
        if isinstance(text_document, file):
            text_document = text_document.read()

        if isinstance(text_document, str):
            text_document = text_document.decode(self.charset, 'ignore')

        # apply filters
        for filter_ in self.filters:
            text_document = filter_(text_document)

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

    def __init__(self, charset='utf-8', min_n=3, max_n=6):
        self.charset = charset
        self.min_n = min_n
        self.max_n = max_n

    def analyze(self, text_document):
        if isinstance(text_document, str):
            text_document = text_document.decode(self.charset, 'ignore')
        text_document = strip_accents(text_document.lower())

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
    """
    Convert a collection of raw documents to a matrix of token counts.

    Parameters
    ----------
    analyzer: WordNGramAnalyzer or CharNGramAnalyzer, optional

    vocabulary: dict, optional
        A dictionary where keys are tokens and values are indices in the matrix.
        This is useful in order to fix the vocabulary in advance.
    """

    def __init__(self, analyzer=DEFAULT_ANALYZER, vocabulary={}):
        self.analyzer = analyzer
        self.vocabulary = vocabulary

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
        matrix = np.zeros((len(docs), len(vocab)), dtype=long)

        for i, doc_dict in enumerate(docs):
            for idx, count in doc_dict.iteritems():
                matrix[i, idx] = count

        return matrix, vocab

    def _build_vectors(self, raw_documents):
        # raw_documents is an iterable so we don't know its length
        vectors = []
        for i, doc in enumerate(raw_documents):
            vector = np.zeros(len(self.vocabulary), dtype=long)
            for token in self.analyzer.analyze(doc):
                try:
                    vector[self.vocabulary[token]] += 1
                except KeyError:
                    # ignore out-of-vocabulary tokens
                    pass
            vectors.append(vector)

        return np.array(vectors)

    def fit(self, raw_documents, y=None):
        _, self.vocabulary = self._build_vectors_and_vocab(raw_documents)

        return self

    def fit_transform(self, raw_documents, y=None):
        """
        Learn the vocabulary dictionary and return the vectors.

        Parameters
        ----------

        raw_documents: iterable
            an iterable which yields either str, unicode or file object

        Returns
        -------
        vectors: array, [n_samples, n_features]
        """
        vectors, self.vocabulary = self._build_vectors_and_vocab(raw_documents)

        return vectors

    def transform(self, raw_documents):
        """
        Return the vectors.

        Parameters
        ----------

        raw_documents: iterable
            an iterable which yields either str, unicode or file object

        Returns
        -------
        vectors: array, [n_samples, n_features]
        """
        if len(self.vocabulary) == 0:
            raise ValueError, "No vocabulary dictionary available..."

        return self._build_vectors(raw_documents)


class TfidfTransformer(BaseEstimator):
    """
    Transform a count matrix to a TF (term-frequency)
    or TF-IDF (term-frequency inverse-document-frequency)
    representation.

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
        """
        Learn the IDF vector (global term weights).

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

    def transform(self, X):
        """
        Transform a count matrix to a TF or TF-IDF representation.

        Parameters
        ----------
        X: array, [n_samples, n_features]
            a matrix of term/token counts

        Returns
        -------
        vectors: array, [n_samples, n_features]
        """
        vectors = X

        if self.use_tf:
            # term-frequencies (normalized counts)
            vectors = vectors / np.sum(X, axis=1).astype(np.float64)[:,np.newaxis]

        if self.use_idf:
            vectors = vectors * self.idf

        return vectors


class Vectorizer(BaseEstimator):
    """
    Convert a collection of raw documents to a matrix.

    Equivalent to CountVectorizer followed by TfidfTransformer.
    """

    def __init__(self, analyzer=DEFAULT_ANALYZER, use_tf=True, use_idf=True):
        self.tc = CountVectorizer(analyzer)
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
            an iterable which yields either str, unicode or file object

        Returns
        -------
        vectors: array, [n_samples, n_features]
        """
        X = self.tc.fit_transform(raw_documents)
        return self.tfidf.fit(X).transform(X)

    def transform(self, raw_documents):
        """
        Return the vectors.

        Parameters
        ----------

        raw_documents: iterable
            an iterable which yields either str, unicode or file object

        Returns
        -------
        vectors: array, [n_samples, n_features]
        """
        X = self.tc.transform(raw_documents)
        return self.tfidf.transform(X)


class HashingVectorizer(object):
    """Compute term frequencies vectors using hashed term space

    See the Hashing-trick related papers referenced by John Langford on this
    page to get a grasp on the usefulness of this representation:

      http://hunch.net/~jl/projects/hash_reps/index.html

    dim is the number of buckets, higher dim means lower collision rate but
    also higher memory requirements and higher processing times on the
    resulting tfidf vectors.

    Documents is a sequence of lists of tokens to initialize the DF estimates.

    TODO handle bigrams in a smart way such as demonstrated here:

      http://streamhacker.com/2010/05/24/text-classification-sentiment-analysis-stopwords-collocations/

    """
    # TODO: implement me using the murmurhash that might be faster: but profile
    # me first :)

    def __init__(self, dim=5000, probes=1, use_idf=True,
                 analyzer=DEFAULT_ANALYZER):
        self.dim = dim
        self.probes = probes
        self.analyzer = analyzer
        self.use_idf = use_idf

        # start counts at one to avoid zero division while
        # computing IDF
        self.df_counts = np.ones(dim, dtype=long)
        self.tf_vectors = None

    def hash_sign(self, token, probe=0):
        """Compute the hash of token with number proble and hashed sign"""
        h = hash(token + (probe * u"#"))
        return abs(h) % self.dim, 1.0 if h % 2 == 0 else -1.0

    def _sample_document(self, text, tf_vector, update_estimates=True):
        """Extract features from text and update running freq estimates"""
        tokens = self.analyzer.analyze(text)
        for token in tokens:
            # TODO add support for cooccurence tokens in a sentence
            # window
            for probe in xrange(self.probes):
                i, incr = self.hash_sign(token, probe)
                tf_vector[i] += incr
        tf_vector /= len(tokens) * self.probes

        if update_estimates and self.use_idf:
            # update the running DF estimate
            self.df_counts += tf_vector != 0.0
        return tf_vector

    def get_idf(self):
        n_samples = float(len(self.tf_vectors))
        return np.log(n_samples / self.df_counts)

    def get_tfidf(self):
        """Compute the TF-log(IDF) vectors of the sampled documents"""
        if self.tf_vectors is None:
            return None
        return self.tf_vectors * self.get_idf()

    def vectorize(self, text_documents):
        """Vectorize a batch of documents in python utf-8 strings or unicode"""
        tf_vectors = np.zeros((len(text_documents), self.dim))
        for i, text in enumerate(text_documents):
            self._sample_document(text, tf_vectors[i])

        if self.tf_vectors is None:
            self.tf_vectors = tf_vectors
        else:
            self.tf_vectors = np.vstack((self.tf_vectors, tf_vectors))

    def vectorize_files(self, document_filepaths):
        """Vectorize a batch of documents stored in utf-8 text files"""
        tf_vectors = np.zeros((len(document_filepaths), self.dim))
        for i, filepath in enumerate(document_filepaths):
            self._sample_document(file(filepath).read(), tf_vectors[i])

        if self.tf_vectors is None:
            self.tf_vectors = tf_vectors
        else:
            self.tf_vectors = np.vstack((self.tf_vectors, tf_vectors))

    def get_vectors(self):
        if self.use_idf:
            return self.get_tfidf()
        else:
            return self.tf_vectors


class SparseHashingVectorizer(object):
    """Compute term freq vectors using hashed term space in a sparse matrix

    The logic is the same as HashingVectorizer but it is possible to use much
    larger dimension vectors without memory issues thanks to the usage of
    scipy.sparse datastructure to store the tf vectors.
    """

    def __init__(self, dim=100000, probes=1, use_idf=True,
                 analyzer=DEFAULT_ANALYZER):
        self.dim = dim
        self.probes = probes
        self.analyzer = analyzer
        self.use_idf = use_idf

        # start counts at one to avoid zero division while
        # computing IDF
        self.df_counts = np.ones(dim, dtype=long)
        self.tf_vectors = None

    def hash_sign(self, token, probe=0):
        h = hash(token + (probe * u"#"))
        return abs(h) % self.dim, 1.0 if h % 2 == 0 else -1.0

    def _sample_document(self, text, tf_vectors, idx=0, update_estimates=True):
        """Extract features from text and update running freq estimates"""

        tokens = self.analyzer.analyze(text)
        counts = defaultdict(lambda: 0.0)
        for token in tokens:
            # TODO add support for cooccurence tokens in a sentence
            # window
            for probe in xrange(self.probes):
                i, incr = self.hash_sign(token, probe)
                counts[i] += incr
        for k, v in counts.iteritems():
            if v == 0.0:
                # can happen if equally frequent conflicting features
                continue
            tf_vectors[idx, k] = v / (len(tokens) * self.probes)

            if update_estimates and self.use_idf:
                # update the running DF estimate
                self.df_counts[k] += 1

    def get_idf(self):
        n_samples = float(self.tf_vectors.shape[0])
        return np.log(n_samples / self.df_counts)

    def get_tfidf(self):
        """Compute the TF-log(IDF) vectors of the sampled documents"""
        coo = self.tf_vectors.tocoo()
        tf_idf = sp.lil_matrix(coo.shape)
        idf = self.get_idf()
        data, row, col = coo.data, coo.row, coo.col
        for i in xrange(len(data)):
            tf_idf[row[i], col[i]] = data[i] * idf[col[i]]
        return tf_idf.tocsr()

    def vectorize(self, text_documents):
        """Vectorize a batch of documents in python utf-8 strings or unicode"""
        tf_vectors = sp.dok_matrix((len(text_documents), self.dim))
        for i, text in enumerate(text_documents):
            self._sample_document(text, tf_vectors, i)

        if self.tf_vectors is None:
            self.tf_vectors = tf_vectors
        else:
            self.tf_vectors = sp.vstack((self.tf_vectors, tf_vectors))

    def vectorize_files(self, document_filepaths):
        """Vectorize a batch of utf-8 text files"""
        tf_vectors = sp.dok_matrix((len(document_filepaths), self.dim))
        for i, filepath in enumerate(document_filepaths):
            self._sample_document(file(filepath).read(), tf_vectors, i)

        if self.tf_vectors is None:
            self.tf_vectors = tf_vectors
        else:
            self.tf_vectors = sp.vstack((self.tf_vectors, tf_vectors))

    def get_vectors(self):
        if self.use_idf:
            return self.get_tfidf()
        else:
            return self.tf_vectors
