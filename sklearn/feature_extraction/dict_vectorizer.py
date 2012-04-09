# Author: Lars Buitinck <L.J.Buitinck@uva.nl>
# License: BSD-style.

from collections import Mapping, Sequence
from operator import itemgetter

import numpy as np
import scipy.sparse as sp

from ..base import BaseEstimator, TransformerMixin
from ..utils import atleast2d_or_csr


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
    names to feature values into Numpy arrays or scipy.sparse matrices for use
    with scikit-learn estimators.

    When feature values are strings, this transformer will do a binary one-hot
    (aka one-of-K) coding: one boolean-valued feature is constructed for each
    of the possible string values that the feature can take on. For instance,
    a feature "f" that can take on the values "ham" and "spam" will become two
    features in the output, one signifying "f=ham", the other "f=spam".

    Features that do not occur in a sample (mapping) will have a zero value
    in the resulting array/matrix.

    Parameters
    ----------
    dtype : callable, optional
        The type of feature values. Passed to Numpy array/scipy.sparse matrix
        constructors as the dtype argument.
    separator: string, optional
        Separator string used when constructing new features for one-hot
        coding.
    sparse: boolean, optional.
        Whether transform should produce scipy.sparse matrices.
        True by default.

    Examples
    --------
    >>> from sklearn.feature_extraction import DictVectorizer
    >>> v = DictVectorizer(sparse=False)
    >>> D = [{'foo': 1, 'bar': 2}, {'foo': 3, 'baz': 1}]
    >>> X = v.fit_transform(D)
    >>> X
    array([[ 2.,  0.,  1.],
           [ 0.,  1.,  3.]])
    >>> v.inverse_transform(X) == \
        [{'bar': 2.0, 'foo': 1.0}, {'baz': 1.0, 'foo': 3.0}]
    True
    >>> v.transform({'foo': 4, 'unseen_feature': 3})
    array([[ 0.,  0.,  4.]])
    """

    def __init__(self, dtype=np.float64, separator="=", sparse=True):
        self.dtype = dtype
        self.separator = separator
        self.sparse = sparse

    def fit(self, X, y=None):
        """Learn a list of feature name -> indices mappings.

        Parameters
        ----------
        X : Mapping or iterable over Mappings
            Dict(s) or Mapping(s) from feature names (arbitrary Python
            objects) to feature values (strings or convertible to dtype).
        y : (ignored)

        Returns
        -------
        self
        """
        X = _tosequence(X)

        # collect all the possible feature names
        feature_names = set()
        for x in X:
            for f, v in x.iteritems():
                if isinstance(v, basestring):
                    f = "%s%s%s" % (f, self.separator, v)
                feature_names.add(f)

        # sort the feature names to define the mapping
        feature_names = sorted(feature_names)
        self.vocabulary_ = dict((f, i) for i, f in enumerate(feature_names))
        self.feature_names_ = feature_names

        return self

    def fit_transform(self, X, y=None):
        """Learn a list of feature name -> indices mappings and transform X.

        Like fit(X) followed by transform(X).

        Parameters
        ----------
        X : Mapping or iterable over Mappings
            Dict(s) or Mapping(s) from feature names (arbitrary Python
            objects) to feature values (strings or convertible to dtype).
        y : (ignored)

        Returns
        -------
        Xa : {array, sparse matrix}
            Feature vectors; always 2-d.
        """
        X = _tosequence(X)
        self.fit(X)
        return self.transform(X)

    def inverse_transform(self, X, dict_type=dict):
        """Transform array or sparse matrix X back to feature mappings.

        X must have been produced by this DictVectorizer's transform or
        fit_transform method; it may only have passed through transformers
        that preserve the number of features and their order.

        In the case of one-hot/one-of-K coding, the constructed feature
        names and values are returned rather than the original ones.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Sample matrix.
        dict_type : callable, optional
            Constructor for feature mappings. Must conform to the
            collections.Mapping API.

        Returns
        -------
        D : list of dict_type objects, length = n_samples
            Feature mappings for the samples in X.
        """
        X = atleast2d_or_csr(X)     # COO matrix is not subscriptable

        names = self.feature_names_
        Xd = [dict_type() for _ in xrange(X.shape[0])]

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
        X : Mapping or iterable over Mappings, length = n_samples
            Dict(s) or Mapping(s) from feature names (arbitrary Python
            objects) to feature values (strings or convertible to dtype).
        y : (ignored)

        Returns
        -------
        Xa : {array, sparse matrix}
            Feature vectors; always 2-d.
        """
        dtype = self.dtype
        vocab = self.vocabulary_

        if self.sparse:
            X = [X] if isinstance(X, Mapping) else X

            i_ind = []
            j_ind = []
            values = []

            for i, x in enumerate(X):
                for f, v in x.iteritems():
                    if isinstance(v, basestring):
                        f = "%s%s%s" % (f, self.separator, v)
                        v = 1
                    try:
                        j = vocab[f]
                        i_ind.append(i)
                        j_ind.append(j)
                        values.append(dtype(v))
                    except KeyError:
                        pass

            shape = (i + 1, len(vocab))
            return sp.coo_matrix((values, (i_ind, j_ind)),
                                 shape=shape, dtype=dtype)

        else:
            X = _tosequence(X)
            Xa = np.zeros((len(X), len(vocab)), dtype=dtype)

            for i, x in enumerate(X):
                for f, v in x.iteritems():
                    if isinstance(v, basestring):
                        f = "%s%s%s" % (f, self.separator, v)
                        v = 1
                    try:
                        Xa[i, vocab[f]] = dtype(v)
                    except KeyError:
                        pass

            return Xa

    def get_feature_names(self):
        """Returns a list of feature names, ordered by their indices.

        If one-of-K coding is applied to categorical features, this will
        include the constructed feature names but not the original ones.
        """
        return self.feature_names_

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

        names = self.feature_names_
        new_vocab = {}
        for i in support:
            new_vocab[names[i]] = len(new_vocab)

        self.vocabulary_ = new_vocab
        self.feature_names_ = [f for f, i in sorted(new_vocab.iteritems(),
                                                    key=itemgetter(1))]

        return self
