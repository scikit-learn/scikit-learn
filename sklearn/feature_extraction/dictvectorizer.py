# Author: Lars Buitinck <L.J.Buitinck@uva.nl>
# License: BSD-style.

from collections import Mapping, Sequence
from operator import itemgetter

import numpy as np
import scipy.sparse as sp

from ..base import BaseEstimator, TransformerMixin


def _tosequence(X):
    """Turn X into a sequence or ndarray, avoiding a copy if possible."""
    if isinstance(X, Mapping):
        return [X]
    elif isinstance(X, (Sequence, np.ndarray)):
        return X
    else:
        return list(X)


class DictVectorizer(BaseEstimator, TransformerMixin):
    """Transforms lists of feature-value mappings to vectors.

    This transformer turns lists of mappings (dict-like objects) of feature
    names to feature values them into Numpy arrays or scipy.sparse matrices
    for use with scikit-learn estimators.

    Features that do not occur in a sample (mapping) will have a zero value
    in the resulting array/matrix.

    Parameters
    ----------
    dtype : callable, optional
        The type of feature values. Passed to Numpy array/scipy.sparse matrix
        constructors as the dtype argument.
    sparse: boolean, optional
        Whether transform should produce scipy.sparse matrices.

    Examples
    --------
    >>> from sklearn.feature_extraction import DictVectorizer
    >>> v = DictVectorizer(sparse=False)
    >>> D = [{'foo': 1, 'bar': 2}, {'foo': 3, 'baz': 1}]
    >>> X = v.fit_transform(D)
    >>> X
    array([[ 1.,  2.,  0.],
           [ 3.,  0.,  1.]])
    >>> v.inverse_transform(X)
    [{'bar': 2.0, 'foo': 1.0}, {'baz': 1.0, 'foo': 3.0}]
    >>> v.transform({'foo': 4, 'quux': 3})
    array([[ 4.,  0.,  0.]])
    """

    def __init__(self, dtype=np.float64, sparse=True):
        self.dtype = dtype
        self.sparse = sparse

    def fit(self, X, y=None):
        """Learn a list of feature name -> indices mappings.

        Parameters
        ----------
        X : Mapping or iterable over Mappings
            Dict(s) or Mapping(s) from feature names (arbitrary Python
            objects) to feature values (must be convertible to dtype).
        y : (ignored)

        Returns
        -------
        self
        """
        X = _tosequence(X)
        vocab = {}

        for x in X:
            for f in x:
                vocab.setdefault(f, len(vocab))

        self.vocabulary_ = vocab

        return self

    def fit_transform(self, X, y=None):
        """Learn a list of feature name -> indices mappings and transform X.

        Like fit(X) followed by transform(X).

        Parameters
        ----------
        X : Mapping or iterable over Mappings
            Dict(s) or Mapping(s) from feature names (arbitrary Python
            objects) to feature values (must be convertible to dtype).
        y : (ignored)

        Returns
        -------
        Xa : {array, sparse matrix}
            Feature vectors; always 2-d.
        """
        X = _tosequence(X)
        self.fit(X)
        return self.transform(X)

    def inverse_transform(self, X, y=None):
        if not sp.issparse(X):
            X = np.atleast_2d(X)

        names = self.get_feature_names()
        Xd = [{} for _ in xrange(X.shape[0])]

        if sp.issparse(X):
            for i, j in zip(*X.nonzero()):
                Xd[i][names[j]] = X[i, j]
        else:
            for i in xrange(X.shape[0]):
                d = Xd[i]
                for j, v in enumerate(X[i, :]):
                    if v != 0:
                        d[names[j]] = X[i, j]

        return Xd

    def transform(self, X, y=None):
        """Transform feature->value dicts to array or sparse matrix.

        Named features not encountered during fit or fit_transform will be
        silently ignored.

        Parameters
        ----------
        X : Mapping or iterable over Mappings
            Dict(s) or Mapping(s) from feature names (arbitrary Python
            objects) to feature values (must be convertible to dtype).
        y : (ignored)

        Returns
        -------
        Xa : {array, sparse matrix}
            Feature vectors; always 2-d.
        """
        X = _tosequence(X)

        zeros = sp.lil_matrix if self.sparse else np.zeros
        dtype = self.dtype
        vocab = self.vocabulary_

        Xa = zeros((len(X), len(vocab)), dtype=dtype)

        for i, x in enumerate(X):
            # The Mapping interface specifies items, but not iteritems;
            # Python 3 will of course solve all our problems.
            try:
                items = x.iteritems()
            except AttributeError:
                items = x.items()

            for f, v in items:
                try:
                    Xa[i, vocab[f]] = dtype(v)
                except KeyError:
                    pass

        return Xa

    def get_feature_names(self):
        """Returns a list of feature names, ordered by their indices."""
        return [f for f, i in sorted(self.vocabulary_.iteritems(),
                                     key=itemgetter(1))]

    def restrict(self, support, indices=False):
        """Restrict the features to those in support.

        Parameters
        ----------
        support : array-like
            Boolean mask or list of indices (as returned by the get_support
            member of feature selectors).
        indices : boolean, optional
            Whether support is a list of indices.
        """
        if not indices:
            support = np.where(support)[0]

        names = self.get_feature_names()
        new_vocab = {}
        for i in support:
            new_vocab[names[i]] = len(new_vocab)

        self.vocabulary_ = new_vocab

        return self
