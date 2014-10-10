# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Mathieu Blondel <mathieu@mblondel.org>
#          Olivier Grisel <olivier.grisel@ensta.org>
#          Andreas Mueller <amueller@ais.uni-bonn.de>
# License: BSD 3 clause

from itertools import chain, combinations
import numbers

import numpy as np
from scipy import sparse

from ..base import BaseEstimator, TransformerMixin
from ..externals import six
from ..utils import check_array
from ..utils import warn_if_not_float
from ..utils.extmath import row_norms
from ..utils.fixes import combinations_with_replacement as combinations_w_r
from ..utils.sparsefuncs_fast import (inplace_csr_row_normalize_l1,
                                      inplace_csr_row_normalize_l2)
from ..utils.sparsefuncs import (inplace_column_scale, mean_variance_axis)

zip = six.moves.zip
map = six.moves.map
range = six.moves.range

__all__ = [
    'Binarizer',
    'KernelCenterer',
    'MinMaxScaler',
    'Normalizer',
    'OneHotEncoder',
    'StandardScaler',
    'add_dummy_feature',
    'binarize',
    'normalize',
    'scale',
]


def _mean_and_std(X, axis=0, with_mean=True, with_std=True):
    """Compute mean and std deviation for centering, scaling.

    Zero valued std components are reset to 1.0 to avoid NaNs when scaling.
    """
    X = np.asarray(X)
    Xr = np.rollaxis(X, axis)

    if with_mean:
        mean_ = Xr.mean(axis=0)
    else:
        mean_ = None

    if with_std:
        std_ = Xr.std(axis=0)
        if isinstance(std_, np.ndarray):
            std_[std_ == 0.0] = 1.0
        elif std_ == 0.:
            std_ = 1.
    else:
        std_ = None

    return mean_, std_


def scale(X, axis=0, with_mean=True, with_std=True, copy=True):
    """Standardize a dataset along any axis

    Center to the mean and component wise scale to unit variance.

    Parameters
    ----------
    X : array-like or CSR matrix.
        The data to center and scale.

    axis : int (0 by default)
        axis used to compute the means and standard deviations along. If 0,
        independently standardize each feature, otherwise (if 1) standardize
        each sample.

    with_mean : boolean, True by default
        If True, center the data before scaling.

    with_std : boolean, True by default
        If True, scale the data to unit variance (or equivalently,
        unit standard deviation).

    copy : boolean, optional, default is True
        set to False to perform inplace row normalization and avoid a
        copy (if the input is already a numpy array or a scipy.sparse
        CSR matrix and if axis is 1).

    Notes
    -----
    This implementation will refuse to center scipy.sparse matrices
    since it would make them non-sparse and would potentially crash the
    program with memory exhaustion problems.

    Instead the caller is expected to either set explicitly
    `with_mean=False` (in that case, only variance scaling will be
    performed on the features of the CSR matrix) or to call `X.toarray()`
    if he/she expects the materialized dense array to fit in memory.

    To avoid memory copy the caller should pass a CSR matrix.

    See also
    --------
    :class:`sklearn.preprocessing.StandardScaler` to perform centering and
    scaling using the ``Transformer`` API (e.g. as part of a preprocessing
    :class:`sklearn.pipeline.Pipeline`)
    """
    if sparse.issparse(X):
        if with_mean:
            raise ValueError(
                "Cannot center sparse matrices: pass `with_mean=False` instead"
                " See docstring for motivation and alternatives.")
        if axis != 0:
            raise ValueError("Can only scale sparse matrix on axis=0, "
                             " got axis=%d" % axis)
        warn_if_not_float(X, estimator='The scale function')
        if not sparse.isspmatrix_csr(X):
            X = X.tocsr()
            copy = False
        if copy:
            X = X.copy()
        _, var = mean_variance_axis(X, axis=0)
        var[var == 0.0] = 1.0
        inplace_column_scale(X, 1 / np.sqrt(var))
    else:
        X = np.asarray(X)
        warn_if_not_float(X, estimator='The scale function')
        mean_, std_ = _mean_and_std(
            X, axis, with_mean=with_mean, with_std=with_std)
        if copy:
            X = X.copy()
        # Xr is a view on the original array that enables easy use of
        # broadcasting on the axis in which we are interested in
        Xr = np.rollaxis(X, axis)
        if with_mean:
            Xr -= mean_
        if with_std:
            Xr /= std_
    return X


class MinMaxScaler(BaseEstimator, TransformerMixin):
    """Standardizes features by scaling each feature to a given range.

    This estimator scales and translates each feature individually such
    that it is in the given range on the training set, i.e. between
    zero and one.

    The standardization is given by::
        X_std = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))
        X_scaled = X_std * (max - min) + min

    where min, max = feature_range.

    This standardization is often used as an alternative to zero mean,
    unit variance scaling.

    Parameters
    ----------
    feature_range: tuple (min, max), default=(0, 1)
        Desired range of transformed data.

    copy : boolean, optional, default is True
        Set to False to perform inplace row normalization and avoid a
        copy (if the input is already a numpy array).

    Attributes
    ----------
    min_ : ndarray, shape (n_features,)
        Per feature adjustment for minimum.

    scale_ : ndarray, shape (n_features,)
        Per feature relative scaling of the data.
    """

    def __init__(self, feature_range=(0, 1), copy=True):
        self.feature_range = feature_range
        self.copy = copy

    def fit(self, X, y=None):
        """Compute the minimum and maximum to be used for later scaling.

        Parameters
        ----------
        X : array-like, shape [n_samples, n_features]
            The data used to compute the per-feature minimum and maximum
            used for later scaling along the features axis.
        """
        X = check_array(X, copy=self.copy, ensure_2d=False)
        warn_if_not_float(X, estimator=self)
        feature_range = self.feature_range
        if feature_range[0] >= feature_range[1]:
            raise ValueError("Minimum of desired feature range must be smaller"
                             " than maximum. Got %s." % str(feature_range))
        data_min = np.min(X, axis=0)
        data_range = np.max(X, axis=0) - data_min
        # Do not scale constant features
        if isinstance(data_range, np.ndarray):
            data_range[data_range == 0.0] = 1.0
        elif data_range == 0.:
            data_range = 1.
        self.scale_ = (feature_range[1] - feature_range[0]) / data_range
        self.min_ = feature_range[0] - data_min * self.scale_
        self.data_range = data_range
        self.data_min = data_min
        return self

    def transform(self, X):
        """Scaling features of X according to feature_range.

        Parameters
        ----------
        X : array-like with shape [n_samples, n_features]
            Input data that will be transformed.
        """
        X = check_array(X, copy=self.copy, ensure_2d=False)
        X *= self.scale_
        X += self.min_
        return X

    def inverse_transform(self, X):
        """Undo the scaling of X according to feature_range.

        Parameters
        ----------
        X : array-like with shape [n_samples, n_features]
            Input data that will be transformed.
        """
        X = check_array(X, copy=self.copy, ensure_2d=False)
        X -= self.min_
        X /= self.scale_
        return X


class StandardScaler(BaseEstimator, TransformerMixin):
    """Standardize features by removing the mean and scaling to unit variance

    Centering and scaling happen independently on each feature by computing
    the relevant statistics on the samples in the training set. Mean and
    standard deviation are then stored to be used on later data using the
    `transform` method.

    Standardization of a dataset is a common requirement for many
    machine learning estimators: they might behave badly if the
    individual feature do not more or less look like standard normally
    distributed data (e.g. Gaussian with 0 mean and unit variance).

    For instance many elements used in the objective function of
    a learning algorithm (such as the RBF kernel of Support Vector
    Machines or the L1 and L2 regularizers of linear models) assume that
    all features are centered around 0 and have variance in the same
    order. If a feature has a variance that is orders of magnitude larger
    that others, it might dominate the objective function and make the
    estimator unable to learn from other features correctly as expected.

    Parameters
    ----------
    with_mean : boolean, True by default
        If True, center the data before scaling.
        This does not work (and will raise an exception) when attempted on
        sparse matrices, because centering them entails building a dense
        matrix which in common use cases is likely to be too large to fit in
        memory.

    with_std : boolean, True by default
        If True, scale the data to unit variance (or equivalently,
        unit standard deviation).

    copy : boolean, optional, default is True
        If False, try to avoid a copy and do inplace scaling instead.
        This is not guaranteed to always work inplace; e.g. if the data is
        not a NumPy array or scipy.sparse CSR matrix, a copy may still be
        returned.

    Attributes
    ----------
    mean_ : array of floats with shape [n_features]
        The mean value for each feature in the training set.

    std_ : array of floats with shape [n_features]
        The standard deviation for each feature in the training set.

    See also
    --------
    :func:`sklearn.preprocessing.scale` to perform centering and
    scaling without using the ``Transformer`` object oriented API

    :class:`sklearn.decomposition.RandomizedPCA` with `whiten=True`
    to further remove the linear correlation across features.
    """

    def __init__(self, copy=True, with_mean=True, with_std=True):
        self.with_mean = with_mean
        self.with_std = with_std
        self.copy = copy

    def fit(self, X, y=None):
        """Compute the mean and std to be used for later scaling.

        Parameters
        ----------
        X : array-like or CSR matrix with shape [n_samples, n_features]
            The data used to compute the mean and standard deviation
            used for later scaling along the features axis.
        """
        X = check_array(X, accept_sparse='csr', copy=self.copy, ensure_2d=False)
        if warn_if_not_float(X, estimator=self):
            X = X.astype(np.float)
        if sparse.issparse(X):
            if self.with_mean:
                raise ValueError(
                    "Cannot center sparse matrices: pass `with_mean=False` "
                    "instead. See docstring for motivation and alternatives.")
            self.mean_ = None

            if self.with_std:
                var = mean_variance_axis(X, axis=0)[1]
                self.std_ = np.sqrt(var)
                self.std_[var == 0.0] = 1.0
            else:
                self.std_ = None
            return self
        else:
            self.mean_, self.std_ = _mean_and_std(
                X, axis=0, with_mean=self.with_mean, with_std=self.with_std)
            return self

    def transform(self, X, y=None, copy=None):
        """Perform standardization by centering and scaling

        Parameters
        ----------
        X : array-like with shape [n_samples, n_features]
            The data used to scale along the features axis.
        """
        copy = copy if copy is not None else self.copy
        X = check_array(X, accept_sparse='csr', copy=copy, ensure_2d=False)
        if warn_if_not_float(X, estimator=self):
            X = X.astype(np.float)
        if sparse.issparse(X):
            if self.with_mean:
                raise ValueError(
                    "Cannot center sparse matrices: pass `with_mean=False` "
                    "instead. See docstring for motivation and alternatives.")
            if self.std_ is not None:
                inplace_column_scale(X, 1 / self.std_)
        else:
            if self.with_mean:
                X -= self.mean_
            if self.with_std:
                X /= self.std_
        return X

    def inverse_transform(self, X, copy=None):
        """Scale back the data to the original representation

        Parameters
        ----------
        X : array-like with shape [n_samples, n_features]
            The data used to scale along the features axis.
        """
        copy = copy if copy is not None else self.copy
        if sparse.issparse(X):
            if self.with_mean:
                raise ValueError(
                    "Cannot uncenter sparse matrices: pass `with_mean=False` "
                    "instead See docstring for motivation and alternatives.")
            if not sparse.isspmatrix_csr(X):
                X = X.tocsr()
                copy = False
            if copy:
                X = X.copy()
            if self.std_ is not None:
                inplace_column_scale(X, self.std_)
        else:
            X = np.asarray(X)
            if copy:
                X = X.copy()
            if self.with_std:
                X *= self.std_
            if self.with_mean:
                X += self.mean_
        return X


class PolynomialFeatures(BaseEstimator, TransformerMixin):
    """Generate polynomial and interaction features.

    Generate a new feature matrix consisting of all polynomial combinations
    of the features with degree less than or equal to the specified degree.
    For example, if an input sample is two dimensional and of the form
    [a, b], the degree-2 polynomial features are [1, a, b, a^2, ab, b^2].

    Parameters
    ----------
    degree : integer
        The degree of the polynomial features. Default = 2.

    interaction_only : boolean, default = False
        If true, only interaction features are produced: features that are
        products of at most ``degree`` *distinct* input features (so not
        ``x[1] ** 2``, ``x[0] * x[2] ** 3``, etc.).

    include_bias : boolean
        If True (default), then include a bias column, the feature in which
        all polynomial powers are zero (i.e. a column of ones - acts as an
        intercept term in a linear model).

    Examples
    --------
    >>> X = np.arange(6).reshape(3, 2)
    >>> X
    array([[0, 1],
           [2, 3],
           [4, 5]])
    >>> poly = PolynomialFeatures(2)
    >>> poly.fit_transform(X)
    array([[ 1,  0,  1,  0,  0,  1],
           [ 1,  2,  3,  4,  6,  9],
           [ 1,  4,  5, 16, 20, 25]])
    >>> poly = PolynomialFeatures(interaction_only=True)
    >>> poly.fit_transform(X)
    array([[ 1,  0,  1,  0],
           [ 1,  2,  3,  6],
           [ 1,  4,  5, 20]])

    Attributes
    ----------

    powers_ :
         powers_[i, j] is the exponent of the jth input in the ith output.

    Notes
    -----
    Be aware that the number of features in the output array scales
    polynomially in the number of features of the input array, and
    exponentially in the degree. High degrees can cause overfitting.

    See :ref:`examples/linear_model/plot_polynomial_interpolation.py
    <example_linear_model_plot_polynomial_interpolation.py>`
    """
    def __init__(self, degree=2, interaction_only=False, include_bias=True):
        self.degree = degree
        self.interaction_only = interaction_only
        self.include_bias = include_bias

    @staticmethod
    def _power_matrix(n_features, degree, interaction_only, include_bias):
        """Compute the matrix of polynomial powers"""
        comb = (combinations if interaction_only else combinations_w_r)
        start = int(not include_bias)
        combn = chain.from_iterable(comb(range(n_features), i)
                                    for i in range(start, degree + 1))
        powers = np.vstack(np.bincount(c, minlength=n_features) for c in combn)
        return powers

    def fit(self, X, y=None):
        """
        Compute the polynomial feature combinations
        """
        n_samples, n_features = check_array(X).shape
        self.powers_ = self._power_matrix(n_features, self.degree,
                                          self.interaction_only,
                                          self.include_bias)
        return self

    def transform(self, X, y=None):
        """Transform data to polynomial features

        Parameters
        ----------
        X : array with shape [n_samples, n_features]
            The data to transform, row by row.

        Returns
        -------
        XP : np.ndarray shape [n_samples, NP]
            The matrix of features, where NP is the number of polynomial
            features generated from the combination of inputs.
        """
        X = check_array(X)
        n_samples, n_features = X.shape

        if n_features != self.powers_.shape[1]:
            raise ValueError("X shape does not match training shape")

        return (X[:, None, :] ** self.powers_).prod(-1)


def normalize(X, norm='l2', axis=1, copy=True):
    """Scale input vectors individually to unit norm (vector length).

    Parameters
    ----------
    X : array or scipy.sparse matrix with shape [n_samples, n_features]
        The data to normalize, element by element.
        scipy.sparse matrices should be in CSR format to avoid an
        un-necessary copy.

    norm : 'l1' or 'l2', optional ('l2' by default)
        The norm to use to normalize each non zero sample (or each non-zero
        feature if axis is 0).

    axis : 0 or 1, optional (1 by default)
        axis used to normalize the data along. If 1, independently normalize
        each sample, otherwise (if 0) normalize each feature.

    copy : boolean, optional, default is True
        set to False to perform inplace row normalization and avoid a
        copy (if the input is already a numpy array or a scipy.sparse
        CSR matrix and if axis is 1).

    See also
    --------
    :class:`sklearn.preprocessing.Normalizer` to perform normalization
    using the ``Transformer`` API (e.g. as part of a preprocessing
    :class:`sklearn.pipeline.Pipeline`)
    """
    if norm not in ('l1', 'l2'):
        raise ValueError("'%s' is not a supported norm" % norm)

    if axis == 0:
        sparse_format = 'csc'
    elif axis == 1:
        sparse_format = 'csr'
    else:
        raise ValueError("'%d' is not a supported axis" % axis)

    X = check_array(X, sparse_format, copy=copy)
    warn_if_not_float(X, 'The normalize function')
    if axis == 0:
        X = X.T

    if sparse.issparse(X):
        if norm == 'l1':
            inplace_csr_row_normalize_l1(X)
        elif norm == 'l2':
            inplace_csr_row_normalize_l2(X)
    else:
        if norm == 'l1':
            norms = np.abs(X).sum(axis=1)
            norms[norms == 0.0] = 1.0
        elif norm == 'l2':
            norms = row_norms(X)
            norms[norms == 0.0] = 1.0
        X /= norms[:, np.newaxis]

    if axis == 0:
        X = X.T

    return X


class Normalizer(BaseEstimator, TransformerMixin):
    """Normalize samples individually to unit norm.

    Each sample (i.e. each row of the data matrix) with at least one
    non zero component is rescaled independently of other samples so
    that its norm (l1 or l2) equals one.

    This transformer is able to work both with dense numpy arrays and
    scipy.sparse matrix (use CSR format if you want to avoid the burden of
    a copy / conversion).

    Scaling inputs to unit norms is a common operation for text
    classification or clustering for instance. For instance the dot
    product of two l2-normalized TF-IDF vectors is the cosine similarity
    of the vectors and is the base similarity metric for the Vector
    Space Model commonly used by the Information Retrieval community.

    Parameters
    ----------
    norm : 'l1' or 'l2', optional ('l2' by default)
        The norm to use to normalize each non zero sample.

    copy : boolean, optional, default is True
        set to False to perform inplace row normalization and avoid a
        copy (if the input is already a numpy array or a scipy.sparse
        CSR matrix).

    Notes
    -----
    This estimator is stateless (besides constructor parameters), the
    fit method does nothing but is useful when used in a pipeline.

    See also
    --------
    :func:`sklearn.preprocessing.normalize` equivalent function
    without the object oriented API
    """

    def __init__(self, norm='l2', copy=True):
        self.norm = norm
        self.copy = copy

    def fit(self, X, y=None):
        """Do nothing and return the estimator unchanged

        This method is just there to implement the usual API and hence
        work in pipelines.
        """
        X = check_array(X, accept_sparse='csr')
        return self

    def transform(self, X, y=None, copy=None):
        """Scale each non zero row of X to unit norm

        Parameters
        ----------
        X : array or scipy.sparse matrix with shape [n_samples, n_features]
            The data to normalize, row by row. scipy.sparse matrices should be
            in CSR format to avoid an un-necessary copy.
        """
        copy = copy if copy is not None else self.copy
        X = check_array(X, accept_sparse='csr')
        return normalize(X, norm=self.norm, axis=1, copy=copy)


def binarize(X, threshold=0.0, copy=True):
    """Boolean thresholding of array-like or scipy.sparse matrix

    Parameters
    ----------
    X : array or scipy.sparse matrix with shape [n_samples, n_features]
        The data to binarize, element by element.
        scipy.sparse matrices should be in CSR or CSC format to avoid an
        un-necessary copy.

    threshold : float, optional (0.0 by default)
        Feature values below or equal to this are replaced by 0, above it by 1.
        Threshold may not be less than 0 for operations on sparse matrices.

    copy : boolean, optional, default is True
        set to False to perform inplace binarization and avoid a copy
        (if the input is already a numpy array or a scipy.sparse CSR / CSC
        matrix and if axis is 1).

    See also
    --------
    :class:`sklearn.preprocessing.Binarizer` to perform binarization
    using the ``Transformer`` API (e.g. as part of a preprocessing
    :class:`sklearn.pipeline.Pipeline`)
    """
    X = check_array(X, accept_sparse=['csr', 'csc'], copy=copy)
    if sparse.issparse(X):
        if threshold < 0:
            raise ValueError('Cannot binarize a sparse matrix with threshold '
                             '< 0')
        cond = X.data > threshold
        not_cond = np.logical_not(cond)
        X.data[cond] = 1
        X.data[not_cond] = 0
        X.eliminate_zeros()
    else:
        cond = X > threshold
        not_cond = np.logical_not(cond)
        X[cond] = 1
        X[not_cond] = 0
    return X


class Binarizer(BaseEstimator, TransformerMixin):
    """Binarize data (set feature values to 0 or 1) according to a threshold

    Values greater than the threshold map to 1, while values less than
    or equal to the threshold map to 0. With the default threshold of 0,
    only positive values map to 1.

    Binarization is a common operation on text count data where the
    analyst can decide to only consider the presence or absence of a
    feature rather than a quantified number of occurrences for instance.

    It can also be used as a pre-processing step for estimators that
    consider boolean random variables (e.g. modelled using the Bernoulli
    distribution in a Bayesian setting).

    Parameters
    ----------
    threshold : float, optional (0.0 by default)
        Feature values below or equal to this are replaced by 0, above it by 1.
        Threshold may not be less than 0 for operations on sparse matrices.

    copy : boolean, optional, default is True
        set to False to perform inplace binarization and avoid a copy (if
        the input is already a numpy array or a scipy.sparse CSR matrix).

    Notes
    -----
    If the input is a sparse matrix, only the non-zero values are subject
    to update by the Binarizer class.

    This estimator is stateless (besides constructor parameters), the
    fit method does nothing but is useful when used in a pipeline.
    """

    def __init__(self, threshold=0.0, copy=True):
        self.threshold = threshold
        self.copy = copy

    def fit(self, X, y=None):
        """Do nothing and return the estimator unchanged

        This method is just there to implement the usual API and hence
        work in pipelines.
        """
        check_array(X, accept_sparse='csr')
        return self

    def transform(self, X, y=None, copy=None):
        """Binarize each element of X

        Parameters
        ----------
        X : array or scipy.sparse matrix with shape [n_samples, n_features]
            The data to binarize, element by element.
            scipy.sparse matrices should be in CSR format to avoid an
            un-necessary copy.
        """
        copy = copy if copy is not None else self.copy
        return binarize(X, threshold=self.threshold, copy=copy)


class KernelCenterer(BaseEstimator, TransformerMixin):
    """Center a kernel matrix

    Let K(x, z) be a kernel defined by phi(x)^T phi(z), where phi is a
    function mapping x to a Hilbert space. KernelCenterer centers (i.e.,
    normalize to have zero mean) the data without explicitly computing phi(x).
    It is equivalent to centering phi(x) with
    sklearn.preprocessing.StandardScaler(with_std=False).
    """

    def fit(self, K, y=None):
        """Fit KernelCenterer

        Parameters
        ----------
        K : numpy array of shape [n_samples, n_samples]
            Kernel matrix.

        Returns
        -------
        self : returns an instance of self.
        """
        K = check_array(K)
        n_samples = K.shape[0]
        self.K_fit_rows_ = np.sum(K, axis=0) / n_samples
        self.K_fit_all_ = self.K_fit_rows_.sum() / n_samples
        return self

    def transform(self, K, y=None, copy=True):
        """Center kernel matrix.

        Parameters
        ----------
        K : numpy array of shape [n_samples1, n_samples2]
            Kernel matrix.

        Returns
        -------
        K_new : numpy array of shape [n_samples1, n_samples2]
        """
        K = check_array(K)
        if copy:
            K = K.copy()

        K_pred_cols = (np.sum(K, axis=1) /
                       self.K_fit_rows_.shape[0])[:, np.newaxis]

        K -= self.K_fit_rows_
        K -= K_pred_cols
        K += self.K_fit_all_

        return K


def add_dummy_feature(X, value=1.0):
    """Augment dataset with an additional dummy feature.

    This is useful for fitting an intercept term with implementations which
    cannot otherwise fit it directly.

    Parameters
    ----------
    X : array or scipy.sparse matrix with shape [n_samples, n_features]
        Data.

    value : float
        Value to use for the dummy feature.

    Returns
    -------

    X : array or scipy.sparse matrix with shape [n_samples, n_features + 1]
        Same data with dummy feature added as first column.

    Examples
    --------

    >>> from sklearn.preprocessing import add_dummy_feature
    >>> add_dummy_feature([[0, 1], [1, 0]])
    array([[ 1.,  0.,  1.],
           [ 1.,  1.,  0.]])
    """
    X = check_array(X, accept_sparse=['csc', 'csr', 'coo'])
    n_samples, n_features = X.shape
    shape = (n_samples, n_features + 1)
    if sparse.issparse(X):
        if sparse.isspmatrix_coo(X):
            # Shift columns to the right.
            col = X.col + 1
            # Column indices of dummy feature are 0 everywhere.
            col = np.concatenate((np.zeros(n_samples), col))
            # Row indices of dummy feature are 0, ..., n_samples-1.
            row = np.concatenate((np.arange(n_samples), X.row))
            # Prepend the dummy feature n_samples times.
            data = np.concatenate((np.ones(n_samples) * value, X.data))
            return sparse.coo_matrix((data, (row, col)), shape)
        elif sparse.isspmatrix_csc(X):
            # Shift index pointers since we need to add n_samples elements.
            indptr = X.indptr + n_samples
            # indptr[0] must be 0.
            indptr = np.concatenate((np.array([0]), indptr))
            # Row indices of dummy feature are 0, ..., n_samples-1.
            indices = np.concatenate((np.arange(n_samples), X.indices))
            # Prepend the dummy feature n_samples times.
            data = np.concatenate((np.ones(n_samples) * value, X.data))
            return sparse.csc_matrix((data, indices, indptr), shape)
        else:
            klass = X.__class__
            return klass(add_dummy_feature(X.tocoo(), value))
    else:
        return np.hstack((np.ones((n_samples, 1)) * value, X))


def _transform_selected(X, transform, selected="all", copy=True):
    """Apply a transform function to portion of selected features

    Parameters
    ----------
    X : array-like or sparse matrix, shape=(n_samples, n_features)
        Dense array or sparse matrix.

    transform : callable
        A callable transform(X) -> X_transformed

    copy : boolean, optional
        Copy X even if it could be avoided.

    selected: "all" or array of indices or mask
        Specify which features to apply the transform to.

    Returns
    -------
    X : array or sparse matrix, shape=(n_samples, n_features_new)
    """
    if selected == "all":
        return transform(X)

    X = check_array(X, accept_sparse='csc', copy=copy)

    if len(selected) == 0:
        return X

    n_features = X.shape[1]
    ind = np.arange(n_features)
    sel = np.zeros(n_features, dtype=bool)
    sel[np.asarray(selected)] = True
    not_sel = np.logical_not(sel)
    n_selected = np.sum(sel)

    if n_selected == 0:
        # No features selected.
        return X
    elif n_selected == n_features:
        # All features selected.
        return transform(X)
    else:
        X_sel = transform(X[:, ind[sel]])
        X_not_sel = X[:, ind[not_sel]]

        if sparse.issparse(X_sel) or sparse.issparse(X_not_sel):
            return sparse.hstack((X_sel, X_not_sel))
        else:
            return np.hstack((X_sel, X_not_sel))


class OneHotEncoder(BaseEstimator, TransformerMixin):
    """Encode categorical integer features using a one-hot aka one-of-K scheme.

    The input to this transformer should be a matrix of integers, denoting
    the values taken on by categorical (discrete) features. The output will be
    a sparse matrix were each column corresponds to one possible value of one
    feature. It is assumed that input features take on values in the range
    [0, n_values).

    This encoding is needed for feeding categorical data to many scikit-learn
    estimators, notably linear models and SVMs with the standard kernels.

    Parameters
    ----------
    n_values : 'auto', int or array of ints
        Number of values per feature.

        - 'auto' : determine value range from training data.
        - int : maximum value for all features.
        - array : maximum value per feature.

    categorical_features: "all" or array of indices or mask
        Specify what features are treated as categorical.

        - 'all' (default): All features are treated as categorical.
        - array of indices: Array of categorical feature indices.
        - mask: Array of length n_features and with dtype=bool.

        Non-categorical features are always stacked to the right of the matrix.

    dtype : number type, default=np.float
        Desired dtype of output.

    sparse : boolean, default=True
        Will return sparse matrix if set True else will return an array.

    Attributes
    ----------
    active_features_ : array
        Indices for active features, meaning values that actually occur
        in the training set. Only available when n_values is ``'auto'``.

    feature_indices_ : array of shape (n_features,)
        Indices to feature ranges.
        Feature ``i`` in the original data is mapped to features
        from ``feature_indices_[i]`` to ``feature_indices_[i+1]``
        (and then potentially masked by `active_features_` afterwards)

    n_values_ : array of shape (n_features,)
        Maximum number of values per feature.

    Examples
    --------
    Given a dataset with three features and two samples, we let the encoder
    find the maximum value per feature and transform the data to a binary
    one-hot encoding.

    >>> from sklearn.preprocessing import OneHotEncoder
    >>> enc = OneHotEncoder()
    >>> enc.fit([[0, 0, 3], [1, 1, 0], [0, 2, 1], \
[1, 0, 2]])  # doctest: +ELLIPSIS
    OneHotEncoder(categorical_features='all', dtype=<... 'float'>,
           n_values='auto', sparse=True)
    >>> enc.n_values_
    array([2, 3, 4])
    >>> enc.feature_indices_
    array([0, 2, 5, 9])
    >>> enc.transform([[0, 1, 1]]).toarray()
    array([[ 1.,  0.,  0.,  1.,  0.,  0.,  1.,  0.,  0.]])

    See also
    --------
    sklearn.feature_extraction.DictVectorizer : performs a one-hot encoding of
      dictionary items (also handles string-valued features).
    sklearn.feature_extraction.FeatureHasher : performs an approximate one-hot
      encoding of dictionary items or strings.
    """
    def __init__(self, n_values="auto", categorical_features="all",
                 dtype=np.float, sparse=True):
        self.n_values = n_values
        self.categorical_features = categorical_features
        self.dtype = dtype
        self.sparse = sparse

    def fit(self, X, y=None):
        """Fit OneHotEncoder to X.

        Parameters
        ----------
        X : array-like, shape=(n_samples, n_feature)
            Input array of type int.

        Returns
        -------
        self
        """
        self.fit_transform(X)
        return self

    def _fit_transform(self, X):
        """Assumes X contains only categorical features."""
        X = check_array(X, dtype=np.int)
        if np.any(X < 0):
            raise ValueError("X needs to contain only non-negative integers.")
        n_samples, n_features = X.shape
        if self.n_values == 'auto':
            n_values = np.max(X, axis=0) + 1
        elif isinstance(self.n_values, numbers.Integral):
            if (np.max(X, axis=0) >= self.n_values).any():
                raise ValueError("Feature out of bounds for n_values=%d"
                                 % self.n_values)
            n_values = np.empty(n_features, dtype=np.int)
            n_values.fill(self.n_values)
        else:
            try:
                n_values = np.asarray(self.n_values, dtype=int)
            except (ValueError, TypeError):
                raise TypeError("Wrong type for parameter `n_values`. Expected"
                                " 'auto', int or array of ints, got %r"
                                % type(X))
            if n_values.ndim < 1 or n_values.shape[0] != X.shape[1]:
                raise ValueError("Shape mismatch: if n_values is an array,"
                                 " it has to be of shape (n_features,).")

        self.n_values_ = n_values
        n_values = np.hstack([[0], n_values])
        indices = np.cumsum(n_values)
        self.feature_indices_ = indices

        column_indices = (X + indices[:-1]).ravel()
        row_indices = np.repeat(np.arange(n_samples, dtype=np.int32),
                                n_features)
        data = np.ones(n_samples * n_features)
        out = sparse.coo_matrix((data, (row_indices, column_indices)),
                                shape=(n_samples, indices[-1]),
                                dtype=self.dtype).tocsr()

        if self.n_values == 'auto':
            mask = np.array(out.sum(axis=0)).ravel() != 0
            active_features = np.where(mask)[0]
            out = out[:, active_features]
            self.active_features_ = active_features

        return out if self.sparse else out.toarray()

    def fit_transform(self, X, y=None):
        """Fit OneHotEncoder to X, then transform X.

        Equivalent to self.fit(X).transform(X), but more convenient and more
        efficient. See fit for the parameters, transform for the return value.
        """
        return _transform_selected(X, self._fit_transform,
                                   self.categorical_features, copy=True)

    def _transform(self, X):
        """Assumes X contains only categorical features."""
        X = check_array(X, dtype=np.int)
        if np.any(X < 0):
            raise ValueError("X needs to contain only non-negative integers.")
        n_samples, n_features = X.shape

        indices = self.feature_indices_
        if n_features != indices.shape[0] - 1:
            raise ValueError("X has different shape than during fitting."
                             " Expected %d, got %d."
                             % (indices.shape[0] - 1, n_features))

        if (np.max(X, axis=0) >= self.n_values_).any():
            raise ValueError("Feature out of bounds. Try setting n_values.")

        column_indices = (X + indices[:-1]).ravel()
        row_indices = np.repeat(np.arange(n_samples, dtype=np.int32),
                                n_features)
        data = np.ones(n_samples * n_features)
        out = sparse.coo_matrix((data, (row_indices, column_indices)),
                                shape=(n_samples, indices[-1]),
                                dtype=self.dtype).tocsr()
        if self.n_values == 'auto':
            out = out[:, self.active_features_]

        return out if self.sparse else out.toarray()

    def transform(self, X):
        """Transform X using one-hot encoding.

        Parameters
        ----------
        X : array-like, shape=(n_samples, n_features)
            Input array of type int.

        Returns
        -------
        X_out : sparse matrix if sparse=True else a 2-d array, dtype=int
            Transformed input.
        """
        return _transform_selected(X, self._transform,
                                   self.categorical_features, copy=True)
