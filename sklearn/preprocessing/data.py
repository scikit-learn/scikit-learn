# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Mathieu Blondel <mathieu@mblondel.org>
#          Olivier Grisel <olivier.grisel@ensta.org>
#          Andreas Mueller <amueller@ais.uni-bonn.de>
#          Eric Martin <eric@ericmart.in>
#          Giorgio Patrini <giorgio.patrini@anu.edu.au>
# License: BSD 3 clause

from __future__ import division

from itertools import chain, combinations
import numbers
import warnings
from itertools import combinations_with_replacement as combinations_w_r

import numpy as np
from scipy import sparse
from scipy import stats

from ..base import BaseEstimator, TransformerMixin
from ..externals import six
from ..externals.six import string_types
from ..utils import check_array
from ..utils.extmath import row_norms
from ..utils.extmath import _incremental_mean_and_var
from ..utils.sparsefuncs_fast import (inplace_csr_row_normalize_l1,
                                      inplace_csr_row_normalize_l2)
from ..utils.sparsefuncs import (inplace_column_scale,
                                 mean_variance_axis, incr_mean_variance_axis,
                                 min_max_axis)
from ..utils.validation import (check_is_fitted, check_random_state,
                                FLOAT_DTYPES)
BOUNDS_THRESHOLD = 1e-7


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
    'QuantileTransformer',
    'add_dummy_feature',
    'binarize',
    'normalize',
    'scale',
    'robust_scale',
    'maxabs_scale',
    'minmax_scale',
    'quantile_transform',
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

    For a comparison of the different scalers, transformers, and normalizers,
    see :ref:`examples/preprocessing/plot_all_scaling.py
    <sphx_glr_auto_examples_preprocessing_plot_all_scaling.py>`.

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

    Examples
    --------
    >>> from sklearn.preprocessing import MinMaxScaler
    >>>
    >>> data = [[-1, 2], [-0.5, 6], [0, 10], [1, 18]]
    >>> scaler = MinMaxScaler()
    >>> print(scaler.fit(data))
    MinMaxScaler(copy=True, feature_range=(0, 1))
    >>> print(scaler.data_max_)
    [  1.  18.]
    >>> print(scaler.transform(data))
    [[ 0.    0.  ]
     [ 0.25  0.25]
     [ 0.5   0.5 ]
     [ 1.    1.  ]]
    >>> print(scaler.transform([[2, 2]]))
    [[ 1.5  0. ]]

    See also
    --------
    minmax_scale: Equivalent function without the estimator API.

    Notes
    -----
    For a comparison of the different scalers, transformers, and normalizers,
    see :ref:`examples/preprocessing/plot_all_scaling.py
    <sphx_glr_auto_examples_preprocessing_plot_all_scaling.py>`.
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
    X : array-like, shape (n_samples, n_features)
        The data.

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

    Notes
    -----
    For a comparison of the different scalers, transformers, and normalizers,
    see :ref:`examples/preprocessing/plot_all_scaling.py
    <sphx_glr_auto_examples_preprocessing_plot_all_scaling.py>`.
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
    copy : boolean, optional, default True
        If False, try to avoid a copy and do inplace scaling instead.
        This is not guaranteed to always work inplace; e.g. if the data is
        not a NumPy array or scipy.sparse CSR matrix, a copy may still be
        returned.

    with_mean : boolean, True by default
        If True, center the data before scaling.
        This does not work (and will raise an exception) when attempted on
        sparse matrices, because centering them entails building a dense
        matrix which in common use cases is likely to be too large to fit in
        memory.

    with_std : boolean, True by default
        If True, scale the data to unit variance (or equivalently,
        unit standard deviation).

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

    Examples
    --------
    >>> from sklearn.preprocessing import StandardScaler
    >>>
    >>> data = [[0, 0], [0, 0], [1, 1], [1, 1]]
    >>> scaler = StandardScaler()
    >>> print(scaler.fit(data))
    StandardScaler(copy=True, with_mean=True, with_std=True)
    >>> print(scaler.mean_)
    [ 0.5  0.5]
    >>> print(scaler.transform(data))
    [[-1. -1.]
     [-1. -1.]
     [ 1.  1.]
     [ 1.  1.]]
    >>> print(scaler.transform([[2, 2]]))
    [[ 3.  3.]]

    See also
    --------
    scale: Equivalent function without the estimator API.

    :class:`sklearn.decomposition.PCA`
        Further removes the linear correlation across features with 'whiten=True'.

    Notes
    -----
    For a comparison of the different scalers, transformers, and normalizers,
    see :ref:`examples/preprocessing/plot_all_scaling.py
    <sphx_glr_auto_examples_preprocessing_plot_all_scaling.py>`.
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

    def transform(self, X, y='deprecated', copy=None):
        """Perform standardization by centering and scaling

        Parameters
        ----------
        X : array-like, shape [n_samples, n_features]
            The data used to scale along the features axis.
        y : (ignored)
            .. deprecated:: 0.19
               This parameter will be removed in 0.21.
        copy : bool, optional (default: None)
            Copy the input X or not.
        """
        if not isinstance(y, string_types) or y != 'deprecated':
            warnings.warn("The parameter y on transform() is "
                          "deprecated since 0.19 and will be removed in 0.21",
                          DeprecationWarning)

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
        copy : bool, optional (default: None)
            Copy the input X or not.

        Returns
        -------
        X_tr : array-like, shape [n_samples, n_features]
            Transformed array.
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
    maxabs_scale: Equivalent function without the estimator API.

    Notes
    -----
    For a comparison of the different scalers, transformers, and normalizers,
    see :ref:`examples/preprocessing/plot_all_scaling.py
    <sphx_glr_auto_examples_preprocessing_plot_all_scaling.py>`.
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

    def transform(self, X):
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
    X : array-like, shape (n_samples, n_features)
        The data.

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

    Notes
    -----
    For a comparison of the different scalers, transformers, and normalizers,
    see :ref:`examples/preprocessing/plot_all_scaling.py
    <sphx_glr_auto_examples_preprocessing_plot_all_scaling.py>`.
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
    sample, depending on the ``axis`` argument) by computing the relevant
    statistics on the samples in the training set. Median and  interquartile
    range are then stored to be used on later data using the ``transform``
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
        This will cause ``transform`` to raise an exception when attempted on
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
    robust_scale: Equivalent function without the estimator API.

    :class:`sklearn.decomposition.PCA`
        Further removes the linear correlation across features with
        'whiten=True'.

    Notes
    -----
    For a comparison of the different scalers, transformers, and normalizers,
    see :ref:`examples/preprocessing/plot_all_scaling.py
    <sphx_glr_auto_examples_preprocessing_plot_all_scaling.py>`.

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

    def transform(self, X):
        """Center and scale the data.

        Can be called on sparse input, provided that ``RobustScaler`` has been
        fitted to dense input and ``with_centering=False``.

        Parameters
        ----------
        X : {array-like, sparse matrix}
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

    For a comparison of the different scalers, transformers, and normalizers,
    see :ref:`examples/preprocessing/plot_all_scaling.py
    <sphx_glr_auto_examples_preprocessing_plot_all_scaling.py>`.

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
        return np.vstack(np.bincount(c, minlength=self.n_input_features_)
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


        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            The data.

        Returns
        -------
        self : instance
        """
        n_samples, n_features = check_array(X).shape
        combinations = self._combinations(n_features, self.degree,
                                          self.interaction_only,
                                          self.include_bias)
        self.n_input_features_ = n_features
        self.n_output_features_ = sum(1 for _ in combinations)
        return self

    def transform(self, X):
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

    Notes
    -----
    For a comparison of the different scalers, transformers, and normalizers,
    see :ref:`examples/preprocessing/plot_all_scaling.py
    <sphx_glr_auto_examples_preprocessing_plot_all_scaling.py>`.

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

    For a comparison of the different scalers, transformers, and normalizers,
    see :ref:`examples/preprocessing/plot_all_scaling.py
    <sphx_glr_auto_examples_preprocessing_plot_all_scaling.py>`.


    See also
    --------
    normalize: Equivalent function without the estimator API.
    """

    def __init__(self, norm='l2', copy=True):
        self.norm = norm
        self.copy = copy

    def fit(self, X, y=None):
        """Do nothing and return the estimator unchanged

        This method is just there to implement the usual API and hence
        work in pipelines.

        Parameters
        ----------
        X : array-like
        """
        X = check_array(X, accept_sparse='csr')
        return self

    def transform(self, X, y='deprecated', copy=None):
        """Scale each non zero row of X to unit norm

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape [n_samples, n_features]
            The data to normalize, row by row. scipy.sparse matrices should be
            in CSR format to avoid an un-necessary copy.
        y : (ignored)
            .. deprecated:: 0.19
               This parameter will be removed in 0.21.
        copy : bool, optional (default: None)
            Copy the input X or not.
        """
        if not isinstance(y, string_types) or y != 'deprecated':
            warnings.warn("The parameter y on transform() is "
                          "deprecated since 0.19 and will be removed in 0.21",
                          DeprecationWarning)

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
    binarize: Equivalent function without the estimator API.
    """

    def __init__(self, threshold=0.0, copy=True):
        self.threshold = threshold
        self.copy = copy

    def fit(self, X, y=None):
        """Do nothing and return the estimator unchanged

        This method is just there to implement the usual API and hence
        work in pipelines.

        Parameters
        ----------
        X : array-like
        """
        check_array(X, accept_sparse='csr')
        return self

    def transform(self, X, y='deprecated', copy=None):
        """Binarize each element of X

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape [n_samples, n_features]
            The data to binarize, element by element.
            scipy.sparse matrices should be in CSR format to avoid an
            un-necessary copy.
        y : (ignored)
            .. deprecated:: 0.19
               This parameter will be removed in 0.21.
        copy : bool
            Copy the input X or not.
        """
        if not isinstance(y, string_types) or y != 'deprecated':
            warnings.warn("The parameter y on transform() is "
                          "deprecated since 0.19 and will be removed in 0.21",
                          DeprecationWarning)

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

    def transform(self, K, y='deprecated', copy=True):
        """Center kernel matrix.

        Parameters
        ----------
        K : numpy array of shape [n_samples1, n_samples2]
            Kernel matrix.
        y : (ignored)
            .. deprecated:: 0.19
               This parameter will be removed in 0.21.
        copy : boolean, optional, default True
            Set to False to perform inplace computation.

        Returns
        -------
        K_new : numpy array of shape [n_samples1, n_samples2]
        """
        if not isinstance(y, string_types) or y != 'deprecated':
            warnings.warn("The parameter y on transform() is "
                          "deprecated since 0.19 and will be removed in 0.21",
                          DeprecationWarning)

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


def _transform_selected(X, transform, selected="all", copy=True):
    """Apply a transform function to portion of selected features

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape [n_samples, n_features]
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
    X = check_array(X, accept_sparse='csc', copy=copy, dtype=FLOAT_DTYPES)

    if isinstance(selected, six.string_types) and selected == "all":
        return transform(X)

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
    a sparse matrix where each column corresponds to one possible value of one
    feature. It is assumed that input features take on values in the range
    [0, n_values).

    This encoding is needed for feeding categorical data to many scikit-learn
    estimators, notably linear models and SVMs with the standard kernels.

    Note: a one-hot encoding of y labels should use a LabelBinarizer
    instead.

    Read more in the :ref:`User Guide <preprocessing_categorical_features>`.

    Parameters
    ----------
    n_values : 'auto', int or array of ints
        Number of values per feature.

        - 'auto' : determine value range from training data.
        - int : number of categorical values per feature.
                Each feature value should be in ``range(n_values)``
        - array : ``n_values[i]`` is the number of categorical values in
                  ``X[:, i]``. Each feature value should be
                  in ``range(n_values[i])``

    categorical_features : "all" or array of indices or mask
        Specify what features are treated as categorical.

        - 'all' (default): All features are treated as categorical.
        - array of indices: Array of categorical feature indices.
        - mask: Array of length n_features and with dtype=bool.

        Non-categorical features are always stacked to the right of the matrix.

    dtype : number type, default=np.float
        Desired dtype of output.

    sparse : boolean, default=True
        Will return sparse matrix if set True else will return an array.

    handle_unknown : str, 'error' or 'ignore'
        Whether to raise an error or ignore if a unknown categorical feature is
        present during transform.

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
    Given a dataset with three features and four samples, we let the encoder
    find the maximum value per feature and transform the data to a binary
    one-hot encoding.

    >>> from sklearn.preprocessing import OneHotEncoder
    >>> enc = OneHotEncoder()
    >>> enc.fit([[0, 0, 3], [1, 1, 0], [0, 2, 1], \
[1, 0, 2]])  # doctest: +ELLIPSIS
    OneHotEncoder(categorical_features='all', dtype=<... 'numpy.float64'>,
           handle_unknown='error', n_values='auto', sparse=True)
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
    sklearn.preprocessing.LabelBinarizer : binarizes labels in a one-vs-all
      fashion.
    sklearn.preprocessing.MultiLabelBinarizer : transforms between iterable of
      iterables and a multilabel format, e.g. a (samples x classes) binary
      matrix indicating the presence of a class label.
    sklearn.preprocessing.LabelEncoder : encodes labels with values between 0
      and n_classes-1.
    """
    def __init__(self, n_values="auto", categorical_features="all",
                 dtype=np.float64, sparse=True, handle_unknown='error'):
        self.n_values = n_values
        self.categorical_features = categorical_features
        self.dtype = dtype
        self.sparse = sparse
        self.handle_unknown = handle_unknown

    def fit(self, X, y=None):
        """Fit OneHotEncoder to X.

        Parameters
        ----------
        X : array-like, shape [n_samples, n_feature]
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
        if (isinstance(self.n_values, six.string_types) and
                self.n_values == 'auto'):
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

        if (isinstance(self.n_values, six.string_types) and
                self.n_values == 'auto'):
            mask = np.array(out.sum(axis=0)).ravel() != 0
            active_features = np.where(mask)[0]
            out = out[:, active_features]
            self.active_features_ = active_features

        return out if self.sparse else out.toarray()

    def fit_transform(self, X, y=None):
        """Fit OneHotEncoder to X, then transform X.

        Equivalent to self.fit(X).transform(X), but more convenient and more
        efficient. See fit for the parameters, transform for the return value.

        Parameters
        ----------
        X : array-like, shape [n_samples, n_feature]
            Input array of type int.
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

        # We use only those categorical features of X that are known using fit.
        # i.e lesser than n_values_ using mask.
        # This means, if self.handle_unknown is "ignore", the row_indices and
        # col_indices corresponding to the unknown categorical feature are
        # ignored.
        mask = (X < self.n_values_).ravel()
        if np.any(~mask):
            if self.handle_unknown not in ['error', 'ignore']:
                raise ValueError("handle_unknown should be either error or "
                                 "unknown got %s" % self.handle_unknown)
            if self.handle_unknown == 'error':
                raise ValueError("unknown categorical feature present %s "
                                 "during transform." % X.ravel()[~mask])

        column_indices = (X + indices[:-1]).ravel()[mask]
        row_indices = np.repeat(np.arange(n_samples, dtype=np.int32),
                                n_features)[mask]
        data = np.ones(np.sum(mask))
        out = sparse.coo_matrix((data, (row_indices, column_indices)),
                                shape=(n_samples, indices[-1]),
                                dtype=self.dtype).tocsr()
        if (isinstance(self.n_values, six.string_types) and
                self.n_values == 'auto'):
            out = out[:, self.active_features_]

        return out if self.sparse else out.toarray()

    def transform(self, X):
        """Transform X using one-hot encoding.

        Parameters
        ----------
        X : array-like, shape [n_samples, n_features]
            Input array of type int.

        Returns
        -------
        X_out : sparse matrix if sparse=True else a 2-d array, dtype=int
            Transformed input.
        """
        return _transform_selected(X, self._transform,
                                   self.categorical_features, copy=True)


class QuantileTransformer(BaseEstimator, TransformerMixin):
    """Transform features using quantiles information.

    This method transforms the features to follow a uniform or a normal
    distribution. Therefore, for a given feature, this transformation tends
    to spread out the most frequent values. It also reduces the impact of
    (marginal) outliers: this is therefore a robust preprocessing scheme.

    The transformation is applied on each feature independently.
    The cumulative density function of a feature is used to project the
    original values. Features values of new/unseen data that fall below
    or above the fitted range will be mapped to the bounds of the output
    distribution. Note that this transform is non-linear. It may distort linear
    correlations between variables measured at the same scale but renders
    variables measured at different scales more directly comparable.

    Read more in the :ref:`User Guide <preprocessing_transformer>`.

    Parameters
    ----------
    n_quantiles : int, optional (default=1000)
        Number of quantiles to be computed. It corresponds to the number
        of landmarks used to discretize the cumulative density function.

    output_distribution : str, optional (default='uniform')
        Marginal distribution for the transformed data. The choices are
        'uniform' (default) or 'normal'.

    ignore_implicit_zeros : bool, optional (default=False)
        Only applies to sparse matrices. If True, the sparse entries of the
        matrix are discarded to compute the quantile statistics. If False,
        these entries are treated as zeros.

    subsample : int, optional (default=1e5)
        Maximum number of samples used to estimate the quantiles for
        computational efficiency. Note that the subsampling procedure may
        differ for value-identical sparse and dense matrices.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by np.random. Note that this is used by subsampling and smoothing
        noise.

    copy : boolean, optional, (default=True)
        Set to False to perform inplace transformation and avoid a copy (if the
        input is already a numpy array).

    Attributes
    ----------
    quantiles_ : ndarray, shape (n_quantiles, n_features)
        The values corresponding the quantiles of reference.

    references_ : ndarray, shape(n_quantiles, )
        Quantiles of references.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.preprocessing import QuantileTransformer
    >>> rng = np.random.RandomState(0)
    >>> X = np.sort(rng.normal(loc=0.5, scale=0.25, size=(25, 1)), axis=0)
    >>> qt = QuantileTransformer(n_quantiles=10, random_state=0)
    >>> qt.fit_transform(X) # doctest: +ELLIPSIS
    array([...])

    See also
    --------
    quantile_transform : Equivalent function without the estimator API.
    StandardScaler : perform standardization that is faster, but less robust
        to outliers.
    RobustScaler : perform robust standardization that removes the influence
        of outliers but does not put outliers and inliers on the same scale.

    Notes
    -----
    For a comparison of the different scalers, transformers, and normalizers,
    see :ref:`examples/preprocessing/plot_all_scaling.py
    <sphx_glr_auto_examples_preprocessing_plot_all_scaling.py>`.
    """

    def __init__(self, n_quantiles=1000, output_distribution='uniform',
                 ignore_implicit_zeros=False, subsample=int(1e5),
                 random_state=None, copy=True):
        self.n_quantiles = n_quantiles
        self.output_distribution = output_distribution
        self.ignore_implicit_zeros = ignore_implicit_zeros
        self.subsample = subsample
        self.random_state = random_state
        self.copy = copy

    def _dense_fit(self, X, random_state):
        """Compute percentiles for dense matrices.

        Parameters
        ----------
        X : ndarray, shape (n_samples, n_features)
            The data used to scale along the features axis.
        """
        if self.ignore_implicit_zeros:
            warnings.warn("'ignore_implicit_zeros' takes effect only with"
                          " sparse matrix. This parameter has no effect.")

        n_samples, n_features = X.shape
        # for compatibility issue with numpy<=1.8.X, references
        # need to be a list scaled between 0 and 100
        references = (self.references_ * 100).tolist()
        self.quantiles_ = []
        for col in X.T:
            if self.subsample < n_samples:
                subsample_idx = random_state.choice(n_samples,
                                                    size=self.subsample,
                                                    replace=False)
                col = col.take(subsample_idx, mode='clip')
            self.quantiles_.append(np.percentile(col, references))
        self.quantiles_ = np.transpose(self.quantiles_)

    def _sparse_fit(self, X, random_state):
        """Compute percentiles for sparse matrices.

        Parameters
        ----------
        X : sparse matrix CSC, shape (n_samples, n_features)
            The data used to scale along the features axis. The sparse matrix
            needs to be nonnegative.
        """
        n_samples, n_features = X.shape

        # for compatibility issue with numpy<=1.8.X, references
        # need to be a list scaled between 0 and 100
        references = list(map(lambda x: x * 100, self.references_))
        self.quantiles_ = []
        for feature_idx in range(n_features):
            column_nnz_data = X.data[X.indptr[feature_idx]:
                                     X.indptr[feature_idx + 1]]
            if len(column_nnz_data) > self.subsample:
                column_subsample = (self.subsample * len(column_nnz_data) //
                                    n_samples)
                if self.ignore_implicit_zeros:
                    column_data = np.zeros(shape=column_subsample,
                                           dtype=X.dtype)
                else:
                    column_data = np.zeros(shape=self.subsample, dtype=X.dtype)
                column_data[:column_subsample] = random_state.choice(
                    column_nnz_data, size=column_subsample, replace=False)
            else:
                if self.ignore_implicit_zeros:
                    column_data = np.zeros(shape=len(column_nnz_data),
                                           dtype=X.dtype)
                else:
                    column_data = np.zeros(shape=n_samples, dtype=X.dtype)
                column_data[:len(column_nnz_data)] = column_nnz_data

            if not column_data.size:
                # if no nnz, an error will be raised for computing the
                # quantiles. Force the quantiles to be zeros.
                self.quantiles_.append([0] * len(references))
            else:
                self.quantiles_.append(
                    np.percentile(column_data, references))
        self.quantiles_ = np.transpose(self.quantiles_)

    def fit(self, X, y=None):
        """Compute the quantiles used for transforming.

        Parameters
        ----------
        X : ndarray or sparse matrix, shape (n_samples, n_features)
            The data used to scale along the features axis. If a sparse
            matrix is provided, it will be converted into a sparse
            ``csc_matrix``. Additionally, the sparse matrix needs to be
            nonnegative if `ignore_implicit_zeros` is False.

        Returns
        -------
        self : object
            Returns self
        """
        if self.n_quantiles <= 0:
            raise ValueError("Invalid value for 'n_quantiles': %d. "
                             "The number of quantiles must be at least one."
                             % self.n_quantiles)

        if self.subsample <= 0:
            raise ValueError("Invalid value for 'subsample': %d. "
                             "The number of subsamples must be at least one."
                             % self.subsample)

        if self.n_quantiles > self.subsample:
            raise ValueError("The number of quantiles cannot be greater than"
                             " the number of samples used. Got {} quantiles"
                             " and {} samples.".format(self.n_quantiles,
                                                       self.subsample))

        X = self._check_inputs(X)
        rng = check_random_state(self.random_state)

        # Create the quantiles of reference
        self.references_ = np.linspace(0, 1, self.n_quantiles,
                                       endpoint=True)
        if sparse.issparse(X):
            self._sparse_fit(X, rng)
        else:
            self._dense_fit(X, rng)

        return self

    def _transform_col(self, X_col, quantiles, inverse):
        """Private function to transform a single feature"""

        if self.output_distribution == 'normal':
            output_distribution = 'norm'
        else:
            output_distribution = self.output_distribution
        output_distribution = getattr(stats, output_distribution)

        # older version of scipy do not handle tuple as fill_value
        # clipping the value before transform solve the issue
        if not inverse:
            lower_bound_x = quantiles[0]
            upper_bound_x = quantiles[-1]
            lower_bound_y = 0
            upper_bound_y = 1
        else:
            lower_bound_x = 0
            upper_bound_x = 1
            lower_bound_y = quantiles[0]
            upper_bound_y = quantiles[-1]
            #  for inverse transform, match a uniform PDF
            X_col = output_distribution.cdf(X_col)
        # find index for lower and higher bounds
        lower_bounds_idx = (X_col - BOUNDS_THRESHOLD <
                            lower_bound_x)
        upper_bounds_idx = (X_col + BOUNDS_THRESHOLD >
                            upper_bound_x)

        if not inverse:
            # Interpolate in one direction and in the other and take the
            # mean. This is in case of repeated values in the features
            # and hence repeated quantiles
            #
            # If we don't do this, only one extreme of the duplicated is
            # used (the upper when we do assending, and the
            # lower for descending). We take the mean of these two
            X_col = .5 * (np.interp(X_col, quantiles, self.references_)
                          - np.interp(-X_col, -quantiles[::-1],
                                      -self.references_[::-1]))
        else:
            X_col = np.interp(X_col, self.references_, quantiles)

        X_col[upper_bounds_idx] = upper_bound_y
        X_col[lower_bounds_idx] = lower_bound_y
        # for forward transform, match the output PDF
        if not inverse:
            X_col = output_distribution.ppf(X_col)
            # find the value to clip the data to avoid mapping to
            # infinity. Clip such that the inverse transform will be
            # consistent
            clip_min = output_distribution.ppf(BOUNDS_THRESHOLD -
                                               np.spacing(1))
            clip_max = output_distribution.ppf(1 - (BOUNDS_THRESHOLD -
                                                    np.spacing(1)))
            X_col = np.clip(X_col, clip_min, clip_max)

        return X_col

    def _check_inputs(self, X, accept_sparse_negative=False):
        """Check inputs before fit and transform"""
        X = check_array(X, accept_sparse='csc', copy=self.copy,
                        dtype=[np.float64, np.float32])
        # we only accept positive sparse matrix when ignore_implicit_zeros is
        # false and that we call fit or transform.
        if (not accept_sparse_negative and not self.ignore_implicit_zeros and
                (sparse.issparse(X) and np.any(X.data < 0))):
            raise ValueError('QuantileTransformer only accepts non-negative'
                             ' sparse matrices.')

        # check the output PDF
        if self.output_distribution not in ('normal', 'uniform'):
            raise ValueError("'output_distribution' has to be either 'normal'"
                             " or 'uniform'. Got '{}' instead.".format(
                                 self.output_distribution))

        return X

    def _check_is_fitted(self, X):
        """Check the inputs before transforming"""
        check_is_fitted(self, 'quantiles_')
        # check that the dimension of X are adequate with the fitted data
        if X.shape[1] != self.quantiles_.shape[1]:
            raise ValueError('X does not have the same number of features as'
                             ' the previously fitted data. Got {} instead of'
                             ' {}.'.format(X.shape[1],
                                           self.quantiles_.shape[1]))

    def _transform(self, X, inverse=False):
        """Forward and inverse transform.

        Parameters
        ----------
        X : ndarray, shape (n_samples, n_features)
            The data used to scale along the features axis.

        inverse : bool, optional (default=False)
            If False, apply forward transform. If True, apply
            inverse transform.

        Returns
        -------
        X : ndarray, shape (n_samples, n_features)
            Projected data
        """

        if sparse.issparse(X):
            for feature_idx in range(X.shape[1]):
                column_slice = slice(X.indptr[feature_idx],
                                     X.indptr[feature_idx + 1])
                X.data[column_slice] = self._transform_col(
                    X.data[column_slice], self.quantiles_[:, feature_idx],
                    inverse)
        else:
            for feature_idx in range(X.shape[1]):
                X[:, feature_idx] = self._transform_col(
                    X[:, feature_idx], self.quantiles_[:, feature_idx],
                    inverse)

        return X

    def transform(self, X):
        """Feature-wise transformation of the data.

        Parameters
        ----------
        X : ndarray or sparse matrix, shape (n_samples, n_features)
            The data used to scale along the features axis. If a sparse
            matrix is provided, it will be converted into a sparse
            ``csc_matrix``. Additionally, the sparse matrix needs to be
            nonnegative if `ignore_implicit_zeros` is False.

        Returns
        -------
        Xt : ndarray or sparse matrix, shape (n_samples, n_features)
            The projected data.
        """
        X = self._check_inputs(X)
        self._check_is_fitted(X)

        return self._transform(X, inverse=False)

    def inverse_transform(self, X):
        """Back-projection to the original space.

        Parameters
        ----------
        X : ndarray or sparse matrix, shape (n_samples, n_features)
            The data used to scale along the features axis. If a sparse
            matrix is provided, it will be converted into a sparse
            ``csc_matrix``. Additionally, the sparse matrix needs to be
            nonnegative if `ignore_implicit_zeros` is False.

        Returns
        -------
        Xt : ndarray or sparse matrix, shape (n_samples, n_features)
            The projected data.
        """
        X = self._check_inputs(X, accept_sparse_negative=True)
        self._check_is_fitted(X)

        return self._transform(X, inverse=True)


def quantile_transform(X, axis=0, n_quantiles=1000,
                       output_distribution='uniform',
                       ignore_implicit_zeros=False,
                       subsample=int(1e5),
                       random_state=None,
                       copy=False):
    """Transform features using quantiles information.

    This method transforms the features to follow a uniform or a normal
    distribution. Therefore, for a given feature, this transformation tends
    to spread out the most frequent values. It also reduces the impact of
    (marginal) outliers: this is therefore a robust preprocessing scheme.

    The transformation is applied on each feature independently.
    The cumulative density function of a feature is used to project the
    original values. Features values of new/unseen data that fall below
    or above the fitted range will be mapped to the bounds of the output
    distribution. Note that this transform is non-linear. It may distort linear
    correlations between variables measured at the same scale but renders
    variables measured at different scales more directly comparable.

    Read more in the :ref:`User Guide <preprocessing_transformer>`.

    Parameters
    ----------
    X : array-like, sparse matrix
        The data to transform.

    axis : int, (default=0)
        Axis used to compute the means and standard deviations along. If 0,
        transform each feature, otherwise (if 1) transform each sample.

    n_quantiles : int, optional (default=1000)
        Number of quantiles to be computed. It corresponds to the number
        of landmarks used to discretize the cumulative density function.

    output_distribution : str, optional (default='uniform')
        Marginal distribution for the transformed data. The choices are
        'uniform' (default) or 'normal'.

    ignore_implicit_zeros : bool, optional (default=False)
        Only applies to sparse matrices. If True, the sparse entries of the
        matrix are discarded to compute the quantile statistics. If False,
        these entries are treated as zeros.

    subsample : int, optional (default=1e5)
        Maximum number of samples used to estimate the quantiles for
        computational efficiency. Note that the subsampling procedure may
        differ for value-identical sparse and dense matrices.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by np.random. Note that this is used by subsampling and smoothing
        noise.

    copy : boolean, optional, (default=True)
        Set to False to perform inplace transformation and avoid a copy (if the
        input is already a numpy array).

    Attributes
    ----------
    quantiles_ : ndarray, shape (n_quantiles, n_features)
        The values corresponding the quantiles of reference.

    references_ : ndarray, shape(n_quantiles, )
        Quantiles of references.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.preprocessing import quantile_transform
    >>> rng = np.random.RandomState(0)
    >>> X = np.sort(rng.normal(loc=0.5, scale=0.25, size=(25, 1)), axis=0)
    >>> quantile_transform(X, n_quantiles=10, random_state=0)
    ... # doctest: +ELLIPSIS
    array([...])

    See also
    --------
    QuantileTransformer : Performs quantile-based scaling using the
        ``Transformer`` API (e.g. as part of a preprocessing
        :class:`sklearn.pipeline.Pipeline`).
    scale : perform standardization that is faster, but less robust
        to outliers.
    robust_scale : perform robust standardization that removes the influence
        of outliers but does not put outliers and inliers on the same scale.

    Notes
    -----
    For a comparison of the different scalers, transformers, and normalizers,
    see :ref:`examples/preprocessing/plot_all_scaling.py
    <sphx_glr_auto_examples_preprocessing_plot_all_scaling.py>`.
    """
    n = QuantileTransformer(n_quantiles=n_quantiles,
                            output_distribution=output_distribution,
                            subsample=subsample,
                            ignore_implicit_zeros=ignore_implicit_zeros,
                            random_state=random_state,
                            copy=copy)
    if axis == 0:
        return n.fit_transform(X)
    elif axis == 1:
        return n.fit_transform(X.T).T
    else:
        raise ValueError("axis should be either equal to 0 or 1. Got"
                         " axis={}".format(axis))
