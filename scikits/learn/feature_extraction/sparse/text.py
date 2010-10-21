
# Authors: Olivier Grisel <olivier.grisel@ensta.org>
#          Mathieu Blondel
#
# License: BSD Style.

from collections import defaultdict
import numpy as np
import scipy.sparse as sp

from ..text import BaseCountVectorizer, BaseTfidfTransformer, BaseVectorizer, \
                   DEFAULT_ANALYZER

class SparseCountVectorizer(BaseCountVectorizer):

    def _init_matrix(self, shape):
        return sp.dok_matrix(shape, dtype=self.dtype)

class SparseTfidfTransformer(BaseTfidfTransformer):

    def fit(self, X, y=None):
        """
        Learn the IDF vector (global term weights).

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
        """
        Transform a count matrix to a TF or TF-IDF representation.

        Parameters
        ----------
        X: sparse matrix, [n_samples, n_features]
            a matrix of term/token counts

        Returns
        -------
        vectors: sparse matrix, [n_samples, n_features]
        """
        X = sp.dok_matrix(X, dtype=np.float64, copy=copy)
        n_samples, n_features = X.shape

        if self.use_tf:
            sums = np.array(X.sum(axis=1)).ravel()
            for doc, token in zip(*X.nonzero()):
                    X[doc, token] /= sums[doc]

        if self.use_idf:
            d = sp.lil_matrix((len(self.idf), len(self.idf)))
            d.setdiag(self.idf)
            X = X * d

        if self.normalize:
            norms = X.multiply(X).sum(axis=1)
            norms = np.sqrt(np.array(norms).ravel())

            for doc, token in zip(*X.nonzero()):
                X[doc, token] /= norms[doc]

        return X

class SparseVectorizer(BaseVectorizer):
    """
    Convert a collection of raw documents to a sparse matrix.

    Equivalent to SparseCountVectorizer followed by SparseTfidfTransformer.
    """

    def __init__(self,
                 analyzer=DEFAULT_ANALYZER,
                 use_tf=True,
                 use_idf=True,
                 normalize=False):
        self.tc = SparseCountVectorizer(analyzer, dtype=np.float64)
        self.tfidf = SparseTfidfTransformer(use_tf, use_idf, normalize)

class SparseHashingVectorizer(object):
    """Compute term freq vectors using hashed term space in a sparse matrix

    The logic is the same as HashingVectorizer but it is possible to use much
    larger dimension vectors without memory issues thanks to the usage of
    scipy.sparse datastructure to store the tf vectors.

    This function requires scipy 0.7 or higher.
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

