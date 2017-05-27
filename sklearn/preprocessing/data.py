# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Mathieu Blondel <mathieu@mblondel.org>
#          Olivier Grisel <olivier.grisel@ensta.org>
#          Andreas Mueller <amueller@ais.uni-bonn.de>
#          Eric Martin <eric@ericmart.in>
#          Giorgio Patrini <giorgio.patrini@anu.edu.au>
# License: BSD 3 clause

from itertools import chain, combinations
import warnings
from itertools import combinations_with_replacement as combinations_w_r

import numpy as np
from scipy import sparse

from ..base import BaseEstimator, TransformerMixin
from ..externals import six
from ..utils import check_array
from ..utils.extmath import row_norms
from ..utils.extmath import _incremental_mean_and_var
from ..utils.fixes import bincount
from ..utils.sparsefuncs_fast import (inplace_csr_row_normalize_l1,
                                      inplace_csr_row_normalize_l2)
from ..utils.sparsefuncs import (inplace_column_scale,
                                 mean_variance_axis, incr_mean_variance_axis,
                                 min_max_axis)
from ..utils.validation import check_is_fitted, FLOAT_DTYPES
from .label import LabelEncoder
from ..utils.fixes import in1d


zip = six.moves.zip
map = six.moves.map
range = six.moves.range

__all__ = [
    'Binarizer',
    'KernelCenterer',
    'MinMaxScaler',
    'MaxAbsScaler',
    'Normalizer',
    'OneHotEncoder',
    'RobustScaler',
    'StandardScaler',
    'add_dummy_feature',
    'binarize',
    'normalize',
    'scale',
    'robust_scale',
    'maxabs_scale',
    'minmax_scale',
]


def _handle_zeros_in_scale(scale, copy=True):
    ''' Makes sure that whenever scale is zero, we handle it correctly.

    This happens in most scalers when we have constant features.'''

    # if we are fitting on 1D arrays, scale might be a scalar
    if np.isscalar(scale):
        if scale == .0:
            scale = 1.
        return scale
    elif isinstance(scale, np.ndarray):
        if copy:
            # New array to avoid side-effects
            scale = scale.copy()
        scale[scale == 0.0] = 1.0
        return scale


def scale(X, axis=0, with_mean=True, with_std=True, copy=True):
    """Standardize a dataset along any axis

    Center to the mean and component wise scale to unit variance.

    Read more in the :ref:`User Guide <preprocessing_scaler>`.

    Parameters
    ----------
    X : {array-like, sparse matrix}
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

    copy : boolean, optional, default True
        set to False to perform inplace row normalization and avoid a
        copy (if the input is already a numpy array or a scipy.sparse
        CSC matrix and if axis is 1).

    Notes
    -----
    This implementation will refuse to center scipy.sparse matrices
    since it would make them non-sparse and would potentially crash the
    program with memory exhaustion problems.

    Instead the caller is expected to either set explicitly
    `with_mean=False` (in that case, only variance scaling will be
    performed on the features of the CSC matrix) or to call `X.toarray()`
    if he/she expects the materialized dense array to fit in memory.

    To avoid memory copy the caller should pass a CSC matrix.

    See also
    --------
    StandardScaler: Performs scaling to unit variance using the``Transformer`` API
        (e.g. as part of a preprocessing :class:`sklearn.pipeline.Pipeline`).
    """  # noqa
    X = check_array(X, accept_sparse='csc', copy=copy, ensure_2d=False,
                    warn_on_dtype=True, estimator='the scale function',
                    dtype=FLOAT_DTYPES)
    if sparse.issparse(X):
        if with_mean:
            raise ValueError(
                "Cannot center sparse matrices: pass `with_mean=False` instead"
                " See docstring for motivation and alternatives.")
        if axis != 0:
            raise ValueError("Can only scale sparse matrix on axis=0, "
                             " got axis=%d" % axis)
        if with_std:
            _, var = mean_variance_axis(X, axis=0)
            var = _handle_zeros_in_scale(var, copy=False)
            inplace_column_scale(X, 1 / np.sqrt(var))
    else:
        X = np.asarray(X)
        if with_mean:
            mean_ = np.mean(X, axis)
        if with_std:
            scale_ = np.std(X, axis)
        # Xr is a view on the original array that enables easy use of
        # broadcasting on the axis in which we are interested in
        Xr = np.rollaxis(X, axis)
        if with_mean:
            Xr -= mean_
            mean_1 = Xr.mean(axis=0)
            # Verify that mean_1 is 'close to zero'. If X contains very
            # large values, mean_1 can also be very large, due to a lack of
            # precision of mean_. In this case, a pre-scaling of the
            # concerned feature is efficient, for instance by its mean or
            # maximum.
            if not np.allclose(mean_1, 0):
                warnings.warn("Numerical issues were encountered "
                              "when centering the data "
                              "and might not be solved. Dataset may "
                              "contain too large values. You may need "
                              "to prescale your features.")
                Xr -= mean_1
        if with_std:
            scale_ = _handle_zeros_in_scale(scale_, copy=False)
            Xr /= scale_
            if with_mean:
                mean_2 = Xr.mean(axis=0)
                # If mean_2 is not 'close to zero', it comes from the fact that
                # scale_ is very small so that mean_2 = mean_1/scale_ > 0, even
                # if mean_1 was close to zero. The problem is thus essentially
                # due to the lack of precision of mean_. A solution is then to
                # subtract the mean again:
                if not np.allclose(mean_2, 0):
                    warnings.warn("Numerical issues were encountered "
                                  "when scaling the data "
                                  "and might not be solved. The standard "
                                  "deviation of the data is probably "
                                  "very close to 0. ")
                    Xr -= mean_2
    return X


class MinMaxScaler(BaseEstimator, TransformerMixin):
    """Transforms features by scaling each feature to a given range.

    This estimator scales and translates each feature individually such
    that it is in the given range on the training set, i.e. between
    zero and one.

    The transformation is given by::

        X_std = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))
        X_scaled = X_std * (max - min) + min

    where min, max = feature_range.

    This transformation is often used as an alternative to zero mean,
    unit variance scaling.

    Read more in the :ref:`User Guide <preprocessing_scaler>`.

    Parameters
    ----------
    feature_range : tuple (min, max), default=(0, 1)
        Desired range of transformed data.

    copy : boolean, optional, default True
        Set to False to perform inplace row normalization and avoid a
        copy (if the input is already a numpy array).

    Attributes
    ----------
    min_ : ndarray, shape (n_features,)
        Per feature adjustment for minimum.

    scale_ : ndarray, shape (n_features,)
        Per feature relative scaling of the data.

        .. versionadded:: 0.17
           *scale_* attribute.

    data_min_ : ndarray, shape (n_features,)
        Per feature minimum seen in the data

        .. versionadded:: 0.17
           *data_min_*

    data_max_ : ndarray, shape (n_features,)
        Per feature maximum seen in the data

        .. versionadded:: 0.17
           *data_max_*

    data_range_ : ndarray, shape (n_features,)
        Per feature range ``(data_max_ - data_min_)`` seen in the data

        .. versionadded:: 0.17
           *data_range_*

    See also
    --------
    minmax_scale: Equivalent function without the object oriented API.
    """

    def __init__(self, feature_range=(0, 1), copy=True):
        self.feature_range = feature_range
        self.copy = copy

    def _reset(self):
        """Reset internal data-dependent state of the scaler, if necessary.

        __init__ parameters are not touched.
        """

        # Checking one attribute is enough, becase they are all set together
        # in partial_fit
        if hasattr(self, 'scale_'):
            del self.scale_
            del self.min_
            del self.n_samples_seen_
            del self.data_min_
            del self.data_max_
            del self.data_range_

    def fit(self, X, y=None):
        """Compute the minimum and maximum to be used for later scaling.

        Parameters
        ----------
        X : array-like, shape [n_samples, n_features]
            The data used to compute the per-feature minimum and maximum
            used for later scaling along the features axis.
        """

        # Reset internal state before fitting
        self._reset()
        return self.partial_fit(X, y)

    def partial_fit(self, X, y=None):
        """Online computation of min and max on X for later scaling.
        All of X is processed as a single batch. This is intended for cases
        when `fit` is not feasible due to very large number of `n_samples`
        or because X is read from a continuous stream.

        Parameters
        ----------
        X : array-like, shape [n_samples, n_features]
            The data used to compute the mean and standard deviation
            used for later scaling along the features axis.

        y : Passthrough for ``Pipeline`` compatibility.
        """
        feature_range = self.feature_range
        if feature_range[0] >= feature_range[1]:
            raise ValueError("Minimum of desired feature range must be smaller"
                             " than maximum. Got %s." % str(feature_range))

        if sparse.issparse(X):
            raise TypeError("MinMaxScaler does no support sparse input. "
                            "You may consider to use MaxAbsScaler instead.")

        X = check_array(X, copy=self.copy, warn_on_dtype=True,
                        estimator=self, dtype=FLOAT_DTYPES)

        data_min = np.min(X, axis=0)
        data_max = np.max(X, axis=0)

        # First pass
        if not hasattr(self, 'n_samples_seen_'):
            self.n_samples_seen_ = X.shape[0]
        # Next steps
        else:
            data_min = np.minimum(self.data_min_, data_min)
            data_max = np.maximum(self.data_max_, data_max)
            self.n_samples_seen_ += X.shape[0]

        data_range = data_max - data_min
        self.scale_ = ((feature_range[1] - feature_range[0]) /
                       _handle_zeros_in_scale(data_range))
        self.min_ = feature_range[0] - data_min * self.scale_
        self.data_min_ = data_min
        self.data_max_ = data_max
        self.data_range_ = data_range
        return self

    def transform(self, X):
        """Scaling features of X according to feature_range.

        Parameters
        ----------
        X : array-like, shape [n_samples, n_features]
            Input data that will be transformed.
        """
        check_is_fitted(self, 'scale_')

        X = check_array(X, copy=self.copy, dtype=FLOAT_DTYPES)

        X *= self.scale_
        X += self.min_
        return X

    def inverse_transform(self, X):
        """Undo the scaling of X according to feature_range.

        Parameters
        ----------
        X : array-like, shape [n_samples, n_features]
            Input data that will be transformed. It cannot be sparse.
        """
        check_is_fitted(self, 'scale_')

        X = check_array(X, copy=self.copy, dtype=FLOAT_DTYPES)

        X -= self.min_
        X /= self.scale_
        return X


def minmax_scale(X, feature_range=(0, 1), axis=0, copy=True):
    """Transforms features by scaling each feature to a given range.

    This estimator scales and translates each feature individually such
    that it is in the given range on the training set, i.e. between
    zero and one.

    The transformation is given by::

        X_std = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))
        X_scaled = X_std * (max - min) + min

    where min, max = feature_range.

    This transformation is often used as an alternative to zero mean,
    unit variance scaling.

    Read more in the :ref:`User Guide <preprocessing_scaler>`.

    .. versionadded:: 0.17
       *minmax_scale* function interface
       to :class:`sklearn.preprocessing.MinMaxScaler`.

    Parameters
    ----------
    feature_range : tuple (min, max), default=(0, 1)
        Desired range of transformed data.

    axis : int (0 by default)
        axis used to scale along. If 0, independently scale each feature,
        otherwise (if 1) scale each sample.

    copy : boolean, optional, default is True
        Set to False to perform inplace scaling and avoid a copy (if the input
        is already a numpy array).

    See also
    --------
    MinMaxScaler: Performs scaling to a given range using the``Transformer`` API
        (e.g. as part of a preprocessing :class:`sklearn.pipeline.Pipeline`).
    """  # noqa
    # Unlike the scaler object, this function allows 1d input.
    # If copy is required, it will be done inside the scaler object.
    X = check_array(X, copy=False, ensure_2d=False, warn_on_dtype=True,
                    dtype=FLOAT_DTYPES)
    original_ndim = X.ndim

    if original_ndim == 1:
        X = X.reshape(X.shape[0], 1)

    s = MinMaxScaler(feature_range=feature_range, copy=copy)
    if axis == 0:
        X = s.fit_transform(X)
    else:
        X = s.fit_transform(X.T).T

    if original_ndim == 1:
        X = X.ravel()

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

    This scaler can also be applied to sparse CSR or CSC matrices by passing
    `with_mean=False` to avoid breaking the sparsity structure of the data.

    Read more in the :ref:`User Guide <preprocessing_scaler>`.

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

    copy : boolean, optional, default True
        If False, try to avoid a copy and do inplace scaling instead.
        This is not guaranteed to always work inplace; e.g. if the data is
        not a NumPy array or scipy.sparse CSR matrix, a copy may still be
        returned.

    Attributes
    ----------
    scale_ : ndarray, shape (n_features,)
        Per feature relative scaling of the data.

        .. versionadded:: 0.17
           *scale_*

    mean_ : array of floats with shape [n_features]
        The mean value for each feature in the training set.

    var_ : array of floats with shape [n_features]
        The variance for each feature in the training set. Used to compute
        `scale_`

    n_samples_seen_ : int
        The number of samples processed by the estimator. Will be reset on
        new calls to fit, but increments across ``partial_fit`` calls.

    See also
    --------
    scale: Equivalent function without the object oriented API.

    :class:`sklearn.decomposition.PCA`
        Further removes the linear correlation across features with 'whiten=True'.
    """  # noqa

    def __init__(self, copy=True, with_mean=True, with_std=True):
        self.with_mean = with_mean
        self.with_std = with_std
        self.copy = copy

    def _reset(self):
        """Reset internal data-dependent state of the scaler, if necessary.

        __init__ parameters are not touched.
        """

        # Checking one attribute is enough, becase they are all set together
        # in partial_fit
        if hasattr(self, 'scale_'):
            del self.scale_
            del self.n_samples_seen_
            del self.mean_
            del self.var_

    def fit(self, X, y=None):
        """Compute the mean and std to be used for later scaling.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape [n_samples, n_features]
            The data used to compute the mean and standard deviation
            used for later scaling along the features axis.

        y : Passthrough for ``Pipeline`` compatibility.
        """

        # Reset internal state before fitting
        self._reset()
        return self.partial_fit(X, y)

    def partial_fit(self, X, y=None):
        """Online computation of mean and std on X for later scaling.
        All of X is processed as a single batch. This is intended for cases
        when `fit` is not feasible due to very large number of `n_samples`
        or because X is read from a continuous stream.

        The algorithm for incremental mean and std is given in Equation 1.5a,b
        in Chan, Tony F., Gene H. Golub, and Randall J. LeVeque. "Algorithms
        for computing the sample variance: Analysis and recommendations."
        The American Statistician 37.3 (1983): 242-247:

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape [n_samples, n_features]
            The data used to compute the mean and standard deviation
            used for later scaling along the features axis.

        y : Passthrough for ``Pipeline`` compatibility.
        """
        X = check_array(X, accept_sparse=('csr', 'csc'), copy=self.copy,
                        warn_on_dtype=True, estimator=self, dtype=FLOAT_DTYPES)

        # Even in the case of `with_mean=False`, we update the mean anyway
        # This is needed for the incremental computation of the var
        # See incr_mean_variance_axis and _incremental_mean_variance_axis

        if sparse.issparse(X):
            if self.with_mean:
                raise ValueError(
                    "Cannot center sparse matrices: pass `with_mean=False` "
                    "instead. See docstring for motivation and alternatives.")
            if self.with_std:
                # First pass
                if not hasattr(self, 'n_samples_seen_'):
                    self.mean_, self.var_ = mean_variance_axis(X, axis=0)
                    self.n_samples_seen_ = X.shape[0]
                # Next passes
                else:
                    self.mean_, self.var_, self.n_samples_seen_ = \
                        incr_mean_variance_axis(X, axis=0,
                                                last_mean=self.mean_,
                                                last_var=self.var_,
                                                last_n=self.n_samples_seen_)
            else:
                self.mean_ = None
                self.var_ = None
        else:
            # First pass
            if not hasattr(self, 'n_samples_seen_'):
                self.mean_ = .0
                self.n_samples_seen_ = 0
                if self.with_std:
                    self.var_ = .0
                else:
                    self.var_ = None

            self.mean_, self.var_, self.n_samples_seen_ = \
                _incremental_mean_and_var(X, self.mean_, self.var_,
                                          self.n_samples_seen_)

        if self.with_std:
            self.scale_ = _handle_zeros_in_scale(np.sqrt(self.var_))
        else:
            self.scale_ = None

        return self

    def transform(self, X, y=None, copy=None):
        """Perform standardization by centering and scaling

        Parameters
        ----------
        X : array-like, shape [n_samples, n_features]
            The data used to scale along the features axis.
        """
        check_is_fitted(self, 'scale_')

        copy = copy if copy is not None else self.copy
        X = check_array(X, accept_sparse='csr', copy=copy, warn_on_dtype=True,
                        estimator=self, dtype=FLOAT_DTYPES)

        if sparse.issparse(X):
            if self.with_mean:
                raise ValueError(
                    "Cannot center sparse matrices: pass `with_mean=False` "
                    "instead. See docstring for motivation and alternatives.")
            if self.scale_ is not None:
                inplace_column_scale(X, 1 / self.scale_)
        else:
            if self.with_mean:
                X -= self.mean_
            if self.with_std:
                X /= self.scale_
        return X

    def inverse_transform(self, X, copy=None):
        """Scale back the data to the original representation

        Parameters
        ----------
        X : array-like, shape [n_samples, n_features]
            The data used to scale along the features axis.
        """
        check_is_fitted(self, 'scale_')

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
            if self.scale_ is not None:
                inplace_column_scale(X, self.scale_)
        else:
            X = np.asarray(X)
            if copy:
                X = X.copy()
            if self.with_std:
                X *= self.scale_
            if self.with_mean:
                X += self.mean_
        return X


class MaxAbsScaler(BaseEstimator, TransformerMixin):
    """Scale each feature by its maximum absolute value.

    This estimator scales and translates each feature individually such
    that the maximal absolute value of each feature in the
    training set will be 1.0. It does not shift/center the data, and
    thus does not destroy any sparsity.

    This scaler can also be applied to sparse CSR or CSC matrices.

    .. versionadded:: 0.17

    Parameters
    ----------
    copy : boolean, optional, default is True
        Set to False to perform inplace scaling and avoid a copy (if the input
        is already a numpy array).

    Attributes
    ----------
    scale_ : ndarray, shape (n_features,)
        Per feature relative scaling of the data.

        .. versionadded:: 0.17
           *scale_* attribute.

    max_abs_ : ndarray, shape (n_features,)
        Per feature maximum absolute value.

    n_samples_seen_ : int
        The number of samples processed by the estimator. Will be reset on
        new calls to fit, but increments across ``partial_fit`` calls.

    See also
    --------
    maxabs_scale: Equivalent function without the object oriented API.
    """

    def __init__(self, copy=True):
        self.copy = copy

    def _reset(self):
        """Reset internal data-dependent state of the scaler, if necessary.

        __init__ parameters are not touched.
        """

        # Checking one attribute is enough, becase they are all set together
        # in partial_fit
        if hasattr(self, 'scale_'):
            del self.scale_
            del self.n_samples_seen_
            del self.max_abs_

    def fit(self, X, y=None):
        """Compute the maximum absolute value to be used for later scaling.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape [n_samples, n_features]
            The data used to compute the per-feature minimum and maximum
            used for later scaling along the features axis.
        """

        # Reset internal state before fitting
        self._reset()
        return self.partial_fit(X, y)

    def partial_fit(self, X, y=None):
        """Online computation of max absolute value of X for later scaling.
        All of X is processed as a single batch. This is intended for cases
        when `fit` is not feasible due to very large number of `n_samples`
        or because X is read from a continuous stream.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape [n_samples, n_features]
            The data used to compute the mean and standard deviation
            used for later scaling along the features axis.

        y : Passthrough for ``Pipeline`` compatibility.
        """
        X = check_array(X, accept_sparse=('csr', 'csc'), copy=self.copy,
                        estimator=self, dtype=FLOAT_DTYPES)

        if sparse.issparse(X):
            mins, maxs = min_max_axis(X, axis=0)
            max_abs = np.maximum(np.abs(mins), np.abs(maxs))
        else:
            max_abs = np.abs(X).max(axis=0)

        # First pass
        if not hasattr(self, 'n_samples_seen_'):
            self.n_samples_seen_ = X.shape[0]
        # Next passes
        else:
            max_abs = np.maximum(self.max_abs_, max_abs)
            self.n_samples_seen_ += X.shape[0]

        self.max_abs_ = max_abs
        self.scale_ = _handle_zeros_in_scale(max_abs)
        return self

    def transform(self, X, y=None):
        """Scale the data

        Parameters
        ----------
        X : {array-like, sparse matrix}
            The data that should be scaled.
        """
        check_is_fitted(self, 'scale_')
        X = check_array(X, accept_sparse=('csr', 'csc'), copy=self.copy,
                        estimator=self, dtype=FLOAT_DTYPES)

        if sparse.issparse(X):
            inplace_column_scale(X, 1.0 / self.scale_)
        else:
            X /= self.scale_
        return X

    def inverse_transform(self, X):
        """Scale back the data to the original representation

        Parameters
        ----------
        X : {array-like, sparse matrix}
            The data that should be transformed back.
        """
        check_is_fitted(self, 'scale_')
        X = check_array(X, accept_sparse=('csr', 'csc'), copy=self.copy,
                        estimator=self, dtype=FLOAT_DTYPES)

        if sparse.issparse(X):
            inplace_column_scale(X, self.scale_)
        else:
            X *= self.scale_
        return X


def maxabs_scale(X, axis=0, copy=True):
    """Scale each feature to the [-1, 1] range without breaking the sparsity.

    This estimator scales each feature individually such
    that the maximal absolute value of each feature in the
    training set will be 1.0.

    This scaler can also be applied to sparse CSR or CSC matrices.

    Parameters
    ----------
    axis : int (0 by default)
        axis used to scale along. If 0, independently scale each feature,
        otherwise (if 1) scale each sample.

    copy : boolean, optional, default is True
        Set to False to perform inplace scaling and avoid a copy (if the input
        is already a numpy array).

    See also
    --------
    MaxAbsScaler: Performs scaling to the [-1, 1] range using the``Transformer`` API
        (e.g. as part of a preprocessing :class:`sklearn.pipeline.Pipeline`).
    """  # noqa
    # Unlike the scaler object, this function allows 1d input.

    # If copy is required, it will be done inside the scaler object.
    X = check_array(X, accept_sparse=('csr', 'csc'), copy=False,
                    ensure_2d=False, dtype=FLOAT_DTYPES)
    original_ndim = X.ndim

    if original_ndim == 1:
        X = X.reshape(X.shape[0], 1)

    s = MaxAbsScaler(copy=copy)
    if axis == 0:
        X = s.fit_transform(X)
    else:
        X = s.fit_transform(X.T).T

    if original_ndim == 1:
        X = X.ravel()

    return X


class RobustScaler(BaseEstimator, TransformerMixin):
    """Scale features using statistics that are robust to outliers.

    This Scaler removes the median and scales the data according to
    the quantile range (defaults to IQR: Interquartile Range).
    The IQR is the range between the 1st quartile (25th quantile)
    and the 3rd quartile (75th quantile).

    Centering and scaling happen independently on each feature (or each
    sample, depending on the `axis` argument) by computing the relevant
    statistics on the samples in the training set. Median and  interquartile
    range are then stored to be used on later data using the `transform`
    method.

    Standardization of a dataset is a common requirement for many
    machine learning estimators. Typically this is done by removing the mean
    and scaling to unit variance. However, outliers can often influence the
    sample mean / variance in a negative way. In such cases, the median and
    the interquartile range often give better results.

    .. versionadded:: 0.17

    Read more in the :ref:`User Guide <preprocessing_scaler>`.

    Parameters
    ----------
    with_centering : boolean, True by default
        If True, center the data before scaling.
        This does not work (and will raise an exception) when attempted on
        sparse matrices, because centering them entails building a dense
        matrix which in common use cases is likely to be too large to fit in
        memory.

    with_scaling : boolean, True by default
        If True, scale the data to interquartile range.

    quantile_range : tuple (q_min, q_max), 0.0 < q_min < q_max < 100.0
        Default: (25.0, 75.0) = (1st quantile, 3rd quantile) = IQR
        Quantile range used to calculate ``scale_``.

        .. versionadded:: 0.18

    copy : boolean, optional, default is True
        If False, try to avoid a copy and do inplace scaling instead.
        This is not guaranteed to always work inplace; e.g. if the data is
        not a NumPy array or scipy.sparse CSR matrix, a copy may still be
        returned.

    Attributes
    ----------
    center_ : array of floats
        The median value for each feature in the training set.

    scale_ : array of floats
        The (scaled) interquartile range for each feature in the training set.

        .. versionadded:: 0.17
           *scale_* attribute.

    See also
    --------
    robust_scale: Equivalent function without the object oriented API.

    :class:`sklearn.decomposition.PCA`
        Further removes the linear correlation across features with
        'whiten=True'.

    Notes
    -----
    See examples/preprocessing/plot_robust_scaling.py for an example.

    https://en.wikipedia.org/wiki/Median_(statistics)
    https://en.wikipedia.org/wiki/Interquartile_range
    """

    def __init__(self, with_centering=True, with_scaling=True,
                 quantile_range=(25.0, 75.0), copy=True):
        self.with_centering = with_centering
        self.with_scaling = with_scaling
        self.quantile_range = quantile_range
        self.copy = copy

    def _check_array(self, X, copy):
        """Makes sure centering is not enabled for sparse matrices."""
        X = check_array(X, accept_sparse=('csr', 'csc'), copy=self.copy,
                        estimator=self, dtype=FLOAT_DTYPES)

        if sparse.issparse(X):
            if self.with_centering:
                raise ValueError(
                    "Cannot center sparse matrices: use `with_centering=False`"
                    " instead. See docstring for motivation and alternatives.")
        return X

    def fit(self, X, y=None):
        """Compute the median and quantiles to be used for scaling.

        Parameters
        ----------
        X : array-like, shape [n_samples, n_features]
            The data used to compute the median and quantiles
            used for later scaling along the features axis.
        """
        if sparse.issparse(X):
            raise TypeError("RobustScaler cannot be fitted on sparse inputs")
        X = self._check_array(X, self.copy)
        if self.with_centering:
            self.center_ = np.median(X, axis=0)

        if self.with_scaling:
            q_min, q_max = self.quantile_range
            if not 0 <= q_min <= q_max <= 100:
                raise ValueError("Invalid quantile range: %s" %
                                 str(self.quantile_range))

            q = np.percentile(X, self.quantile_range, axis=0)
            self.scale_ = (q[1] - q[0])
            self.scale_ = _handle_zeros_in_scale(self.scale_, copy=False)
        return self

    def transform(self, X, y=None):
        """Center and scale the data

        Parameters
        ----------
        X : array-like
            The data used to scale along the specified axis.
        """
        if self.with_centering:
            check_is_fitted(self, 'center_')
        if self.with_scaling:
            check_is_fitted(self, 'scale_')
        X = self._check_array(X, self.copy)

        if sparse.issparse(X):
            if self.with_scaling:
                inplace_column_scale(X, 1.0 / self.scale_)
        else:
            if self.with_centering:
                X -= self.center_
            if self.with_scaling:
                X /= self.scale_
        return X

    def inverse_transform(self, X):
        """Scale back the data to the original representation

        Parameters
        ----------
        X : array-like
            The data used to scale along the specified axis.
        """
        if self.with_centering:
            check_is_fitted(self, 'center_')
        if self.with_scaling:
            check_is_fitted(self, 'scale_')
        X = self._check_array(X, self.copy)

        if sparse.issparse(X):
            if self.with_scaling:
                inplace_column_scale(X, self.scale_)
        else:
            if self.with_scaling:
                X *= self.scale_
            if self.with_centering:
                X += self.center_
        return X


def robust_scale(X, axis=0, with_centering=True, with_scaling=True,
                 quantile_range=(25.0, 75.0), copy=True):
    """Standardize a dataset along any axis

    Center to the median and component wise scale
    according to the interquartile range.

    Read more in the :ref:`User Guide <preprocessing_scaler>`.

    Parameters
    ----------
    X : array-like
        The data to center and scale.

    axis : int (0 by default)
        axis used to compute the medians and IQR along. If 0,
        independently scale each feature, otherwise (if 1) scale
        each sample.

    with_centering : boolean, True by default
        If True, center the data before scaling.

    with_scaling : boolean, True by default
        If True, scale the data to unit variance (or equivalently,
        unit standard deviation).

    quantile_range : tuple (q_min, q_max), 0.0 < q_min < q_max < 100.0
        Default: (25.0, 75.0) = (1st quantile, 3rd quantile) = IQR
        Quantile range used to calculate ``scale_``.

        .. versionadded:: 0.18

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
    `with_centering=False` (in that case, only variance scaling will be
    performed on the features of the CSR matrix) or to call `X.toarray()`
    if he/she expects the materialized dense array to fit in memory.

    To avoid memory copy the caller should pass a CSR matrix.

    See also
    --------
    RobustScaler: Performs centering and scaling using the ``Transformer`` API
        (e.g. as part of a preprocessing :class:`sklearn.pipeline.Pipeline`).
    """
    s = RobustScaler(with_centering=with_centering, with_scaling=with_scaling,
                     quantile_range=quantile_range, copy=copy)
    if axis == 0:
        return s.fit_transform(X)
    else:
        return s.fit_transform(X.T).T


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
    array([[  1.,   0.,   1.,   0.,   0.,   1.],
           [  1.,   2.,   3.,   4.,   6.,   9.],
           [  1.,   4.,   5.,  16.,  20.,  25.]])
    >>> poly = PolynomialFeatures(interaction_only=True)
    >>> poly.fit_transform(X)
    array([[  1.,   0.,   1.,   0.],
           [  1.,   2.,   3.,   6.],
           [  1.,   4.,   5.,  20.]])

    Attributes
    ----------
    powers_ : array, shape (n_output_features, n_input_features)
        powers_[i, j] is the exponent of the jth input in the ith output.

    n_input_features_ : int
        The total number of input features.

    n_output_features_ : int
        The total number of polynomial output features. The number of output
        features is computed by iterating over all suitably sized combinations
        of input features.

    Notes
    -----
    Be aware that the number of features in the output array scales
    polynomially in the number of features of the input array, and
    exponentially in the degree. High degrees can cause overfitting.

    See :ref:`examples/linear_model/plot_polynomial_interpolation.py
    <sphx_glr_auto_examples_linear_model_plot_polynomial_interpolation.py>`
    """
    def __init__(self, degree=2, interaction_only=False, include_bias=True):
        self.degree = degree
        self.interaction_only = interaction_only
        self.include_bias = include_bias

    @staticmethod
    def _combinations(n_features, degree, interaction_only, include_bias):
        comb = (combinations if interaction_only else combinations_w_r)
        start = int(not include_bias)
        return chain.from_iterable(comb(range(n_features), i)
                                   for i in range(start, degree + 1))

    @property
    def powers_(self):
        check_is_fitted(self, 'n_input_features_')

        combinations = self._combinations(self.n_input_features_, self.degree,
                                          self.interaction_only,
                                          self.include_bias)
        return np.vstack(bincount(c, minlength=self.n_input_features_)
                         for c in combinations)

    def get_feature_names(self, input_features=None):
        """
        Return feature names for output features

        Parameters
        ----------
        input_features : list of string, length n_features, optional
            String names for input features if available. By default,
            "x0", "x1", ... "xn_features" is used.

        Returns
        -------
        output_feature_names : list of string, length n_output_features

        """
        powers = self.powers_
        if input_features is None:
            input_features = ['x%d' % i for i in range(powers.shape[1])]
        feature_names = []
        for row in powers:
            inds = np.where(row)[0]
            if len(inds):
                name = " ".join("%s^%d" % (input_features[ind], exp)
                                if exp != 1 else input_features[ind]
                                for ind, exp in zip(inds, row[inds]))
            else:
                name = "1"
            feature_names.append(name)
        return feature_names

    def fit(self, X, y=None):
        """
        Compute number of output features.
        """
        n_samples, n_features = check_array(X).shape
        combinations = self._combinations(n_features, self.degree,
                                          self.interaction_only,
                                          self.include_bias)
        self.n_input_features_ = n_features
        self.n_output_features_ = sum(1 for _ in combinations)
        return self

    def transform(self, X, y=None):
        """Transform data to polynomial features

        Parameters
        ----------
        X : array-like, shape [n_samples, n_features]
            The data to transform, row by row.

        Returns
        -------
        XP : np.ndarray shape [n_samples, NP]
            The matrix of features, where NP is the number of polynomial
            features generated from the combination of inputs.
        """
        check_is_fitted(self, ['n_input_features_', 'n_output_features_'])

        X = check_array(X, dtype=FLOAT_DTYPES)
        n_samples, n_features = X.shape

        if n_features != self.n_input_features_:
            raise ValueError("X shape does not match training shape")

        # allocate output data
        XP = np.empty((n_samples, self.n_output_features_), dtype=X.dtype)

        combinations = self._combinations(n_features, self.degree,
                                          self.interaction_only,
                                          self.include_bias)
        for i, c in enumerate(combinations):
            XP[:, i] = X[:, c].prod(1)

        return XP


def normalize(X, norm='l2', axis=1, copy=True, return_norm=False):
    """Scale input vectors individually to unit norm (vector length).

    Read more in the :ref:`User Guide <preprocessing_normalization>`.

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape [n_samples, n_features]
        The data to normalize, element by element.
        scipy.sparse matrices should be in CSR format to avoid an
        un-necessary copy.

    norm : 'l1', 'l2', or 'max', optional ('l2' by default)
        The norm to use to normalize each non zero sample (or each non-zero
        feature if axis is 0).

    axis : 0 or 1, optional (1 by default)
        axis used to normalize the data along. If 1, independently normalize
        each sample, otherwise (if 0) normalize each feature.

    copy : boolean, optional, default True
        set to False to perform inplace row normalization and avoid a
        copy (if the input is already a numpy array or a scipy.sparse
        CSR matrix and if axis is 1).

    return_norm : boolean, default False
        whether to return the computed norms

    Returns
    -------
    X : {array-like, sparse matrix}, shape [n_samples, n_features]
        Normalized input X.

    norms : array, shape [n_samples] if axis=1 else [n_features]
        An array of norms along given axis for X.
        When X is sparse, a NotImplementedError will be raised
        for norm 'l1' or 'l2'.

    See also
    --------
    Normalizer: Performs normalization using the ``Transformer`` API
        (e.g. as part of a preprocessing :class:`sklearn.pipeline.Pipeline`).
    """
    if norm not in ('l1', 'l2', 'max'):
        raise ValueError("'%s' is not a supported norm" % norm)

    if axis == 0:
        sparse_format = 'csc'
    elif axis == 1:
        sparse_format = 'csr'
    else:
        raise ValueError("'%d' is not a supported axis" % axis)

    X = check_array(X, sparse_format, copy=copy,
                    estimator='the normalize function', dtype=FLOAT_DTYPES)
    if axis == 0:
        X = X.T

    if sparse.issparse(X):
        if return_norm and norm in ('l1', 'l2'):
            raise NotImplementedError("return_norm=True is not implemented "
                                      "for sparse matrices with norm 'l1' "
                                      "or norm 'l2'")
        if norm == 'l1':
            inplace_csr_row_normalize_l1(X)
        elif norm == 'l2':
            inplace_csr_row_normalize_l2(X)
        elif norm == 'max':
            _, norms = min_max_axis(X, 1)
            norms_elementwise = norms.repeat(np.diff(X.indptr))
            mask = norms_elementwise != 0
            X.data[mask] /= norms_elementwise[mask]
    else:
        if norm == 'l1':
            norms = np.abs(X).sum(axis=1)
        elif norm == 'l2':
            norms = row_norms(X)
        elif norm == 'max':
            norms = np.max(X, axis=1)
        norms = _handle_zeros_in_scale(norms, copy=False)
        X /= norms[:, np.newaxis]

    if axis == 0:
        X = X.T

    if return_norm:
        return X, norms
    else:
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

    Read more in the :ref:`User Guide <preprocessing_normalization>`.

    Parameters
    ----------
    norm : 'l1', 'l2', or 'max', optional ('l2' by default)
        The norm to use to normalize each non zero sample.

    copy : boolean, optional, default True
        set to False to perform inplace row normalization and avoid a
        copy (if the input is already a numpy array or a scipy.sparse
        CSR matrix).

    Notes
    -----
    This estimator is stateless (besides constructor parameters), the
    fit method does nothing but is useful when used in a pipeline.

    See also
    --------
    normalize: Equivalent function without the object oriented API.
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
        X : {array-like, sparse matrix}, shape [n_samples, n_features]
            The data to normalize, row by row. scipy.sparse matrices should be
            in CSR format to avoid an un-necessary copy.
        """
        copy = copy if copy is not None else self.copy
        X = check_array(X, accept_sparse='csr')
        return normalize(X, norm=self.norm, axis=1, copy=copy)


def binarize(X, threshold=0.0, copy=True):
    """Boolean thresholding of array-like or scipy.sparse matrix

    Read more in the :ref:`User Guide <preprocessing_binarization>`.

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape [n_samples, n_features]
        The data to binarize, element by element.
        scipy.sparse matrices should be in CSR or CSC format to avoid an
        un-necessary copy.

    threshold : float, optional (0.0 by default)
        Feature values below or equal to this are replaced by 0, above it by 1.
        Threshold may not be less than 0 for operations on sparse matrices.

    copy : boolean, optional, default True
        set to False to perform inplace binarization and avoid a copy
        (if the input is already a numpy array or a scipy.sparse CSR / CSC
        matrix and if axis is 1).

    See also
    --------
    Binarizer: Performs binarization using the ``Transformer`` API
        (e.g. as part of a preprocessing :class:`sklearn.pipeline.Pipeline`).
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

    Read more in the :ref:`User Guide <preprocessing_binarization>`.

    Parameters
    ----------
    threshold : float, optional (0.0 by default)
        Feature values below or equal to this are replaced by 0, above it by 1.
        Threshold may not be less than 0 for operations on sparse matrices.

    copy : boolean, optional, default True
        set to False to perform inplace binarization and avoid a copy (if
        the input is already a numpy array or a scipy.sparse CSR matrix).

    Notes
    -----
    If the input is a sparse matrix, only the non-zero values are subject
    to update by the Binarizer class.

    This estimator is stateless (besides constructor parameters), the
    fit method does nothing but is useful when used in a pipeline.

    See also
    --------
    binarize: Equivalent function without the object oriented API.
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
        X : {array-like, sparse matrix}, shape [n_samples, n_features]
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

    Read more in the :ref:`User Guide <kernel_centering>`.
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
        K = check_array(K, dtype=FLOAT_DTYPES)
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

        copy : boolean, optional, default True
            Set to False to perform inplace computation.

        Returns
        -------
        K_new : numpy array of shape [n_samples1, n_samples2]
        """
        check_is_fitted(self, 'K_fit_all_')

        K = check_array(K, copy=copy, dtype=FLOAT_DTYPES)

        K_pred_cols = (np.sum(K, axis=1) /
                       self.K_fit_rows_.shape[0])[:, np.newaxis]

        K -= self.K_fit_rows_
        K -= K_pred_cols
        K += self.K_fit_all_

        return K

    @property
    def _pairwise(self):
        return True


def add_dummy_feature(X, value=1.0):
    """Augment dataset with an additional dummy feature.

    This is useful for fitting an intercept term with implementations which
    cannot otherwise fit it directly.

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape [n_samples, n_features]
        Data.

    value : float
        Value to use for the dummy feature.

    Returns
    -------

    X : {array, sparse matrix}, shape [n_samples, n_features + 1]
        Same data with dummy feature added as first column.

    Examples
    --------

    >>> from sklearn.preprocessing import add_dummy_feature
    >>> add_dummy_feature([[0, 1], [1, 0]])
    array([[ 1.,  0.,  1.],
           [ 1.,  1.,  0.]])
    """
    X = check_array(X, accept_sparse=['csc', 'csr', 'coo'], dtype=FLOAT_DTYPES)
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


def _apply_selected(X, transform, selected="all", dtype=np.float, copy=True,
                    return_val=True):
    """Apply a function to portion of selected features

    Parameters
    ----------
    X : {array, sparse matrix}, shape [n_samples, n_features]
        Dense array or sparse matrix.
    transform : callable
        A callable transform(X) -> X_transformed
    dtype : dtype
        Cast outputs to this data type
    copy : boolean, optional
        Copy X even if it could be avoided.
    selected: "all" or array of indices or mask
        Specify which features to apply the transform to.
    return_val : boolean, optional
        Whether to return the transformed matrix. If not set `None` is
        returned.

    Returns
    -------
    X : array or sparse matrix, shape=(n_samples, n_features_new)
    """
    if copy:
        X = X.copy()

    if isinstance(selected, six.string_types) and selected == "all":
        X_trans = transform(X)
        return X_trans.astype(dtype) if return_val else None

    if len(selected) == 0:
        return X.astype(dtype) if return_val else None

    n_features = X.shape[1]
    sel = np.zeros(n_features, dtype=bool)
    sel[np.asarray(selected)] = True
    not_sel = np.logical_not(sel)
    n_selected = np.sum(sel)

    if n_selected == 0:
        # No features selected.
        return X.astype(dtype) if return_val else None
    elif n_selected == n_features:
        # All features selected.
        X_trans = transform(X)
        return X_trans.astype(dtype) if return_val else None
    else:
        ind = np.arange(n_features)
        X_sel = transform(X[:, ind[sel]])

        if return_val:
            X_sel = X_sel.astype(dtype)
            X_not_sel = X[:, ind[not_sel]].astype(dtype)
            if sparse.issparse(X_sel) or sparse.issparse(X_not_sel):
                return sparse.hstack((X_sel, X_not_sel))
            else:
                return np.hstack((X_sel, X_not_sel))


class OneHotEncoder(BaseEstimator, TransformerMixin):
    """Encode categorical features using a one-hot aka one-of-K scheme.

    The input to this transformer should be a matrix of integers or strings,
    denoting the values taken on by categorical (discrete) features. The
    output will be a sparse matrix where each column corresponds to one
    possible value of one feature.

    This encoding is needed for feeding categorical data to many scikit-learn
    estimators, notably linear models and SVMs with the standard kernels.

    Note: a one-hot encoding of y labels should use a LabelBinarizer
    instead.

    Read more in the :ref:`User Guide <preprocessing_categorical_features>`.

    Parameters
    ----------
    values : 'auto' or List[List[objects]]
        - 'auto' (default) : Encoded values are those found in training data.
        - list of lists : values for feature ``i`` are in ``values[i]``

    categorical_features : "all" or array of indices or mask
        Specify what features are treated as categorical.

        - 'all' (default): All features are treated as categorical.
        - array of indices: Array of categorical feature indices.
        - mask: Array of length n_features and with dtype=bool.

        Non-categorical features are always stacked to the right of the matrix.

    dtype : number type, default=np.float64
        Desired dtype of output.

    sparse : boolean, default=True
        Will return sparse matrix if set True else will return an array.

    handle_unknown : {'error', 'error-strict', 'ignore'}
        - 'ignore': Ignore all unknown feature values.
        - 'error': Raise an error when the value of an integer feature is more
            than the maximum value seen during fit or less than zero, or when
            the value of a non-integer feature was unseen during ``fit``.
        - 'error-strict': Raise an error when the value of a feature is unseen
            during ``fit``.

    Attributes
    ----------
    feature_index_range_ : array, shape (n_feature, 2)
        ``feature_index_range_[i]`` specifies the range of column indices
        occupied by the input feature `i` in the one-hot encoded array.

    one_hot_feature_index_ : array, shape (n_features_new,)
        ``one_hot_feature_index_[i]`` specifies which feature of the input
        is encoded by column ``i`` in the one-hot encoded array.

    categories_ : array, shape (n_features_new,)
        np.object array containing the category encoded in each feature
        of the output (or None for non-categorical features)

    n_values_ : array of shape (n_features,)
        Number of encoded categories per feature. Has value `0` for
        non-categorical features.

    Examples
    --------
    Given a dataset with two features and three samples, we let the encoder
    find the categories in each feature and transform the data to a binary
    one-hot encoding.

    >>> from sklearn.preprocessing import OneHotEncoder
    >>> enc = OneHotEncoder()
    >>> enc.fit(np.array([['cat', 4], ['mouse', 15], ['dog', 17]], dtype='O'))\
        # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
        OneHotEncoder(categorical_features='all',
           dtype=<... 'numpy.float64'>, handle_unknown='error', n_values=None,
           sparse=True, values='auto')
    >>> enc.n_values_
    array([3, 3])
    >>> enc.feature_index_range_
    array([[0, 3],
           [3, 6]])
    >>> enc.one_hot_feature_index_
    array([0, 0, 0, 1, 1, 1])
    >>> (enc.categories_ ==
    ...  np.array(['cat', 'dog', 'mouse', 4, 15, 17], dtype='O')).all()
    True
    >>> enc.transform([['dog', 4]]).toarray()
    array([[ 0.,  1.,  0.,  1.,  0.,  0.]])

    See also
    --------
    sklearn.feature_extraction.DictVectorizer : performs a one-hot encoding of
      dictionary items (also handles string-valued features).
    sklearn.feature_extraction.FeatureHasher : performs an approximate one-hot
      encoding of dictionary items or strings.
    sklearn.preprocessing.LabelBinarizer : binarizes labels in a one-vs-all
      fashion.
    sklearn.preprocessing.MultiLabelBinarizer : transforms between iterable of
      iterables and a multilabel format, e.g. a (samples x classes) binary
      matrix indicating the presence of a class label.
    sklearn.preprocessing.LabelEncoder : encodes labels with values between 0
      and n_classes-1.
    """
    def __init__(self, values='auto', categorical_features="all",
                 n_values=None, dtype=np.float64, sparse=True,
                 handle_unknown='error'):
        self.values = values
        self.categorical_features = categorical_features
        self.dtype = dtype
        self.sparse = sparse
        self.handle_unknown = handle_unknown
        self.n_values = n_values

    def fit(self, X, y=None):
        """Fit the OneHotEncoder to X.

        Parameters
        ----------
        X : array-like, shape [n_samples, n_feature]
            Array of ints or strings or both.

        Returns
        -------
        self
        """
        if self.handle_unknown not in ['ignore', 'error', 'error-strict']:
            template = ("handle_unknown should be either 'error', "
                        "'error-strict', or 'ignore', got %s")
            raise ValueError(template % self.handle_unknown)
        elif self.handle_unknown == 'error':
            warnings.warn('The behavior of handle_unknown="error" is '
                          'deprecated and will be changed to be the same '
                          'as "error-strict" in version 0.21', FutureWarning)

        X = check_array(X, dtype=None, accept_sparse='csc', copy=False)
        n_samples, n_features = X.shape
        self.n_features_ = n_features

        _apply_selected(X, self._fit, dtype=self.dtype, return_val=False,
                        selected=self.categorical_features, copy=False)

        # Record which columns of output data
        # correspond to each column of input data
        self.feature_index_range_ = np.zeros((n_features, 2), dtype=np.int)

        if isinstance(self.categorical_features, six.string_types) and \
           self.categorical_features == "all":
            categorical = np.ones(n_features, dtype=bool)
        else:
            categorical = np.zeros(n_features, dtype=bool)
            categorical[np.asarray(self.categorical_features)] = True

        start, end = 0, 0
        for i_cat, i_feat in enumerate(np.where(categorical)[0]):
            if np.isscalar(self._values) and self.handle_unknown == 'error':
                end = start + self._n_active_features_[i_cat]
            else:
                end = start + len(self._label_encoders[i_cat].classes_)
            self.feature_index_range_[i_feat] = start, end
            start = end
        num_cat_cols = np.sum(categorical)
        non_cat_indices = np.arange(start, start + n_features - num_cat_cols)
        self.feature_index_range_[~categorical, 0] = non_cat_indices
        self.feature_index_range_[~categorical, 1] = non_cat_indices + 1

        # Record which column of input data corresponds
        # to each column of output data
        n_cats = np.diff(self.feature_index_range_, axis=1).ravel()
        inp_order = np.argsort(self.feature_index_range_[:, 0])
        self.one_hot_feature_index_ = np.repeat(inp_order, n_cats[inp_order])

        # Count categories per feature
        self.n_values_ = n_cats.copy()
        self.n_values_[~categorical] = 0

        # Store categories for each output feature
        if num_cat_cols == 0:
            cats = []
        else:
            cats = np.concatenate([le.classes_ for le in self._label_encoders])
        if hasattr(self, '_active_features_'):
            cats = cats[self._active_features_]
        self.categories_ = np.hstack([cats, len(non_cat_indices) * [None]])

        return self

    def _check_values(self, values, n_features):
        """Verify that the input `values` is valid

        Raises ValueError or TypeError for bad `values`.
        Assume that integers or lists of integers have been
        converted to lists of arrays before getting here.
        This should run after `_initialize_values`.
        """
        error_msg = ("`values` should be 'auto', an integer, "
                     "a list of integers or a list of list")
        if isinstance(values, six.string_types):
            if values != 'auto':
                raise ValueError(error_msg)
        elif isinstance(values, list) or isinstance(values, np.ndarray):
            if len(values) != n_features:
                raise ValueError("Shape mismatch: if values is a list,"
                                 " it has to be of length (n_features).")

            # All entries must be either arrays or lists here
            if any([np.isscalar(val) for val in values]):
                raise ValueError(error_msg)
        else:
            raise TypeError(error_msg)

    def _initialize_values(self):
        """Standardize the `values` input

        Output is either a string or a list of arrays.
        """
        if self.n_values is not None:
            warnings.warn('`n_values` has been renamed to `values`.'
                          'The parameter `n_values` has been deprecated '
                          'and will be removed in version 0.21; use the '
                          'parameter `values` instead and specify the '
                          'expected values for each feature.', FutureWarning)
            values = self.n_values
        else:
            values = self.values

        # Convert `int` and `Sequence[int]` inputs to `List[Array[int]]`
        if (not isinstance(values, six.string_types) and
                np.isscalar(values)):
            warnings.warn('Integer input to `values` is deprecated and'
                          ' will be removed in version 0.21. Specify a '
                          'list of allowed values for each feature instead.',
                          FutureWarning)
            values = np.ones(self.n_features_cat_, dtype=int) * values
        if (not isinstance(values, six.string_types) and
                np.isscalar(values[0])):
            warnings.warn('List of integer input to `values` is deprecated and'
                          ' will be removed in version 0.21. Specify a '
                          'list of allowed values for each feature instead.',
                          FutureWarning)
            values = [np.arange(v, dtype=np.int) for v in values]

        return values

    def _fit(self, X):
        """Assumes `X` contains only categorical features"""
        n_samples, n_features = X.shape
        self.n_features_cat_ = n_features
        self._label_encoders = [LabelEncoder() for i in range(n_features)]

        self._values = self._initialize_values()
        self._check_values(self._values, n_features)

        # Fit on categorical features in the data
        _auto_int_classes = n_features * [None]
        for i in range(n_features):
            le = self._label_encoders[i]

            if np.isscalar(self._values) and self.handle_unknown == 'error':
                # For integer features, allow integers between
                # 0 and column max. The transform will still only
                # return dummy columns for integers present in training data.
                if (not isinstance(X[0, i], six.string_types) and
                        int(X[0, i]) == X[0, i]):
                    _auto_int_classes[i] = np.unique(X[:, i]).astype(int)
                    if np.min(_auto_int_classes[i]) < 0:
                        msg = ('Column %s has value(s) less than zero; all '
                               'integer columns must have minimum value '
                               '0 when value="auto" and '
                               'handle_unknown="error".')
                        raise ValueError(msg)
                    n_classes = np.max(_auto_int_classes[i]) + 1
                    le.fit(np.arange(n_classes))
                else:
                    le.fit(X[:, i])
            elif np.isscalar(self._values):
                le.fit(X[:, i])
            else:
                le.fit(self._values[i])

        if np.isscalar(self._values) and self.handle_unknown == 'error':
            # Record which integer features were present in training
            # data so we can restrict output columns.
            active_features = []
            for i_col, int_classes in enumerate(_auto_int_classes):
                if int_classes is None:
                    n_classes = len(self._label_encoders[i_col].classes_)
                    active_features.append(np.ones(n_classes, dtype=bool))
                else:
                    n_classes = max(self._label_encoders[i_col].classes_) + 1
                    this_col_mask = np.zeros(n_classes, dtype=bool)
                    this_col_mask[int_classes] = True
                    active_features.append(this_col_mask)
            self._n_active_features_ = np.array([a.sum()
                                                 for a in active_features])
            self._active_features_ = np.where(np.hstack(active_features))[0]

    def transform(self, X, y=None):
        """Encode the selected categorical features using the one-hot scheme.

        Parameters
        ----------
        X : array-like, shape [n_samples, n_feature]
            Array of ints or strings or both.

        Returns
        -------
        out : array, shape[n_samples, n_features_new]
            `X` encoded using the one-hot scheme. Will be a CSR sparse
            array if `self.sparse` is True.
        """
        X = check_array(X, accept_sparse='csc', dtype=None, copy=False)
        if X.shape[1] != self.n_features_:
            raise ValueError("Input data must have %s "
                             "features." % self.n_features_)

        return _apply_selected(X, self._transform, dtype=self.dtype,
                               selected=self.categorical_features, copy=False)

    def _transform(self, X):
        """Assumes `X` contains only categorical features."""
        n_samples, n_features = X.shape
        X_int = np.zeros_like(X, dtype=np.int32)

        # Recode all columns of input data as integers
        if self.handle_unknown in ['error', 'error-strict']:
            for i, le in enumerate(self._label_encoders):
                try:
                    X_int[:, i] = le.transform(X[:, i])
                except ValueError as err:
                    orig_msg = str(err)
                    if not orig_msg.startswith('y contains'):
                        raise
                    else:
                        msg = 'Column %d %s' % (i, orig_msg[2:])
                    raise ValueError(msg)
            mask = slice(None)
        else:
            X_mask = np.ones_like(X, dtype=np.bool)
            for i, le in enumerate(self._label_encoders):
                valid_mask = in1d(X[:, i], le.classes_)
                if not np.all(valid_mask):
                    X_mask[:, i] = valid_mask
                    X_int[valid_mask, i] = le.transform(X[valid_mask, i])
                else:
                    X_int[:, i] = le.transform(X[:, i])
            mask = X_mask.ravel()

        # Convert integer columns to sparse array of binary indicators
        n_values = [0] + [le.classes_.shape[0] for le in self._label_encoders]
        indices = np.cumsum(n_values)

        column_indices = (X_int + indices[:-1]).ravel()[mask]
        row_indices = np.repeat(np.arange(n_samples, dtype=np.int32),
                                n_features)[mask]
        data = np.ones(len(row_indices), dtype=self.dtype)

        out = sparse.coo_matrix((data, (row_indices, column_indices)),
                                shape=(n_samples, indices[-1]),
                                dtype=self.dtype).tocsr()

        if np.isscalar(self._values) and self.handle_unknown == 'error':
            out = out[:, self._active_features_]

        return out if self.sparse else out.toarray()

    @property
    def active_features_(self):
        warnings.warn('The property `active_features_` is deprecated and'
                      ' will be removed in version 0.21', FutureWarning)
        if not hasattr(self, '_active_features_'):
            raise AttributeError("'OneHotEncoder' object has no attribute "
                                 "'active_features_'.")
        return self._active_features_

    @property
    def feature_indices_(self):
        # This is very similar to the current attribute
        # `feature_index_range_`, but only applies to the
        # subset of categorical features.
        warnings.warn('The property `feature_indices_` is deprecated and'
                      ' will be removed in version 0.21', FutureWarning)
        if not hasattr(self, '_label_encoders'):
            raise AttributeError("'OneHotEncoder' object has no attribute "
                                 "'feature_indices_'.")
        n_categories = [len(le.classes_) for le in self._label_encoders]
        return np.cumsum([0] + n_categories)
