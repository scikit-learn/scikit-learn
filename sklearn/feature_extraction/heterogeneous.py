from scipy import sparse
import numpy as np

from ..base import BaseEstimator, TransformerMixin
from ..pipeline import _fit_one_transformer, _fit_transform_one, _transform_one
from ..externals.joblib import Parallel, delayed
from ..externals.six import iteritems, string_types


def _getitem(X, column):
    """
    Get feature column from input data (array, dataframe, dict)

    """
    # check whether we have column names or integers
    if isinstance(column, list):
        if all(isinstance(col, int) for col in column):
            column_names = False
        elif all(isinstance(col, string_types) for col in column):
            column_names = True
        else:
            raise ValueError("no valid 'column' type")
    elif isinstance(column, int):
        column_names = False
    elif isinstance(column, string_types):
        column_names = True
    else:
        raise ValueError("no valid 'column' type")

    if column_names:
        return X[column]
    else:
        if hasattr(X, 'iloc'):
            # pandas dataframes
            return X.iloc[:, column]
        else:
            # numpy arrays
            return X[:, column]


class ColumnTransformer(BaseEstimator, TransformerMixin):
    """Applies transformers to columns of a array / dataframe / dict.

    This estimator applies transformer objects to columns or fields of the
    input, then concatenates the results. This is useful for heterogeneous or
    columnar data, to combine several feature extraction mechanisms into a
    single transformer.

    Read more in the :ref:`User Guide <column_transformer>`.

    Parameters
    ----------
    transformers : list of tuples
        List of (name, transformer, column) tuples specifying the transformer
        objects to be applied to subsets of the data. The columns can be
        specified as a scalar or list (for multiple columns) of integer or
        string values. Integers are interpreted as the positional columns,
        strings as the keys of `X`.

    n_jobs : int, optional
        Number of jobs to run in parallel (default 1).

    transformer_weights : dict, optional
        Multiplicative weights for features per transformer.
        Keys are transformer names, values the weights.

    Examples
    --------
    >>> from sklearn.feature_extraction import ColumnTransformer
    >>> from sklearn.preprocessing import Normalizer
    >>> union = ColumnTransformer(
    ...     [("norm1", Normalizer(norm='l1'), 'subset1'),
    ...      ("norm2", Normalizer(norm='l1'), 'subset2')])
    >>> X = {'subset1': [[0., 1.], [2., 2.]], 'subset2': [[1., 1.], [0., 1.]]}
    >>> union.fit_transform(X)    # doctest: +NORMALIZE_WHITESPACE
    array([[ 0. ,  1. ,  0.5,  0.5],
           [ 0.5,  0.5,  0. ,  1. ]])

    """
    def __init__(self, transformers, n_jobs=1, transformer_weights=None):
        self.transformers = transformers
        self.n_jobs = n_jobs
        self.transformer_weights = transformer_weights

    def _iter(self):
        """Generate (name, trans, column, weight) tuples
        """
        get_weight = (self.transformer_weights or {}).get
        return ((name, trans, column, get_weight(name))
                for name, trans, column in self.transformers)

    def get_feature_names(self):
        """Get feature names from all transformers.

        Returns
        -------
        feature_names : list of strings
            Names of the features produced by transform.
        """
        feature_names = []
        for name, trans, column in self.transformers:
            if not hasattr(trans, 'get_feature_names'):
                raise AttributeError("Transformer %s does not provide"
                                     " get_feature_names." % str(name))
            feature_names.extend([name + "__" + f for f in
                                  trans.get_feature_names()])
        return feature_names

    def get_params(self, deep=True):
        if not deep:
            return super(ColumnTransformer, self).get_params(deep=False)
        else:
            out = dict()
            for name, trans, _ in self.transformers:
                for key, value in iteritems(trans.get_params(deep=True)):
                    out['%s__%s' % (name, key)] = value
            out.update(super(ColumnTransformer, self).get_params(deep=False))
            return out

    def fit(self, X, y=None):
        """Fit all transformers using X.

        Parameters
        ----------
        X : array-like or sparse matrix, shape (n_samples, n_features)
            Input data, used to fit transformers.
        """
        transformers = Parallel(n_jobs=self.n_jobs)(
            delayed(_fit_one_transformer)(trans, _getitem(X, column), y)
            for name, trans, column in self.transformers)
        self._update_transformers(transformers)
        return self

    def fit_transform(self, X, y=None, **fit_params):
        """Fit all transformers using X, transform the data and concatenate
        results.

        Parameters
        ----------
        X : array-like or sparse matrix, shape (n_samples, n_features)
            Input data to be transformed.

        Returns
        -------
        X_t : array-like or sparse matrix, shape (n_samples, sum_n_components)
            hstack of results of transformers. sum_n_components is the
            sum of n_components (output dimension) over transformers.
        """
        result = Parallel(n_jobs=self.n_jobs)(
            delayed(_fit_transform_one)(trans, weight, _getitem(X, column), y,
                                        **fit_params)
            for name, trans, column, weight in self._iter())

        Xs, transformers = zip(*result)
        self._update_transformers(transformers)
        if any(sparse.issparse(f) for f in Xs):
            Xs = sparse.hstack(Xs).tocsr()
        else:
            Xs = np.hstack(Xs)
        return Xs

    def transform(self, X):
        """Transform X separately by each transformer, concatenate results.

        Parameters
        ----------
        X : array-like or sparse matrix, shape (n_samples, n_features)
            Input data to be transformed.

        Returns
        -------
        X_t : array-like or sparse matrix, shape (n_samples, sum_n_components)
            hstack of results of transformers. sum_n_components is the
            sum of n_components (output dimension) over transformers.
        """
        Xs = Parallel(n_jobs=self.n_jobs)(
            delayed(_transform_one)(trans, weight, _getitem(X, column))
            for name, trans, column, weight in self._iter())
        if any(sparse.issparse(f) for f in Xs):
            Xs = sparse.hstack(Xs).tocsr()
        else:
            Xs = np.hstack(Xs)
        return Xs

    def _update_transformers(self, transformers):
        # use a dict constructor instead of a dict comprehension for python2.6
        self.transformers = [
            (name, new, column) for ((name, old, column), new)
            in zip(self.transformers, transformers)]
