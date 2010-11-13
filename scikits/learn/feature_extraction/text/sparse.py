
# Authors: Olivier Grisel <olivier.grisel@ensta.org>
#          Mathieu Blondel
#
# License: BSD Style.

import numpy as np
import scipy.sparse as sp

from .dense import BaseCountVectorizer
from .dense import BaseTfidfTransformer
from .dense import BaseVectorizer
from .dense import DEFAULT_ANALYZER

from ...preprocessing.sparse import Normalizer


class CountVectorizer(BaseCountVectorizer):
    """Convert a collection of raw documents to a matrix of token counts

    This implementation produces a sparse representation of the counts using
    scipy.sparse.coo_matrix.

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

    def _term_count_dicts_to_matrix(self, term_count_dicts, vocabulary):

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


class TfidfTransformer(BaseTfidfTransformer):

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


class Vectorizer(BaseVectorizer):
    """Convert a collection of raw documents to a sparse matrix

    Equivalent to CountVectorizer followed by TfidfTransformer.
    """

    def __init__(self, analyzer=DEFAULT_ANALYZER, max_df=1.0,
                 max_features=None, use_tf=True, use_idf=True):
        self.tc = CountVectorizer(analyzer, max_df=max_df,
                                  max_features=max_features,
                                  dtype=np.float64)
        self.tfidf = TfidfTransformer(use_tf, use_idf)


