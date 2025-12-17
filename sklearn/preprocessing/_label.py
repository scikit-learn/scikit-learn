# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import array
import itertools
import warnings
from collections import defaultdict
from numbers import Integral

import numpy as np
import scipy.sparse as sp

from sklearn.base import BaseEstimator, TransformerMixin, _fit_context
from sklearn.utils import column_or_1d
from sklearn.utils._array_api import (
    _convert_to_numpy,
    _find_matching_floating_dtype,
    _is_numpy_namespace,
    _isin,
    device,
    get_namespace,
    get_namespace_and_device,
    indexing_dtype,
    xpx,
)
from sklearn.utils._encode import _encode, _unique
from sklearn.utils._param_validation import Interval, validate_params
from sklearn.utils.multiclass import type_of_target, unique_labels
from sklearn.utils.sparsefuncs import min_max_axis
from sklearn.utils.validation import _num_samples, check_array, check_is_fitted

__all__ = [
    "LabelBinarizer",
    "LabelEncoder",
    "MultiLabelBinarizer",
    "label_binarize",
]


class LabelEncoder(TransformerMixin, BaseEstimator, auto_wrap_output_keys=None):
    """Encode target labels with value between 0 and n_classes-1.

    This transformer should be used to encode target values, *i.e.* `y`, and
    not the input `X`.

    Read more in the :ref:`User Guide <preprocessing_targets>`.

    .. versionadded:: 0.12

    Attributes
    ----------
    classes_ : ndarray of shape (n_classes,)
        Holds the label for each class.

    See Also
    --------
    OrdinalEncoder : Encode categorical features using an ordinal encoding
        scheme.
    OneHotEncoder : Encode categorical features as a one-hot numeric array.

    Examples
    --------
    `LabelEncoder` can be used to normalize labels.

    >>> from sklearn.preprocessing import LabelEncoder
    >>> le = LabelEncoder()
    >>> le.fit([1, 2, 2, 6])
    LabelEncoder()
    >>> le.classes_
    array([1, 2, 6])
    >>> le.transform([1, 1, 2, 6])
    array([0, 0, 1, 2]...)
    >>> le.inverse_transform([0, 0, 1, 2])
    array([1, 1, 2, 6])

    It can also be used to transform non-numerical labels (as long as they are
    hashable and comparable) to numerical labels.

    >>> le = LabelEncoder()
    >>> le.fit(["paris", "paris", "tokyo", "amsterdam"])
    LabelEncoder()
    >>> list(le.classes_)
    [np.str_('amsterdam'), np.str_('paris'), np.str_('tokyo')]
    >>> le.transform(["tokyo", "tokyo", "paris"])
    array([2, 2, 1]...)
    >>> list(le.inverse_transform([2, 2, 1]))
    [np.str_('tokyo'), np.str_('tokyo'), np.str_('paris')]
    """

    def fit(self, y):
        """Fit label encoder.

        Parameters
        ----------
        y : array-like of shape (n_samples,)
            Target values.

        Returns
        -------
        self : returns an instance of self.
            Fitted label encoder.
        """
        y = column_or_1d(y, warn=True)
        self.classes_ = _unique(y)
        return self

    def fit_transform(self, y):
        """Fit label encoder and return encoded labels.

        Parameters
        ----------
        y : array-like of shape (n_samples,)
            Target values.

        Returns
        -------
        y : array-like of shape (n_samples,)
            Encoded labels.
        """
        y = column_or_1d(y, warn=True)
        self.classes_, y = _unique(y, return_inverse=True)
        return y

    def transform(self, y):
        """Transform labels to normalized encoding.

        Parameters
        ----------
        y : array-like of shape (n_samples,)
            Target values.

        Returns
        -------
        y : array-like of shape (n_samples,)
            Labels as normalized encodings.
        """
        check_is_fitted(self)
        xp, _ = get_namespace(y)
        y = column_or_1d(y, dtype=self.classes_.dtype, warn=True)
        # transform of empty array is empty array
        if _num_samples(y) == 0:
            return xp.asarray([])

        return _encode(y, uniques=self.classes_)

    def inverse_transform(self, y):
        """Transform labels back to original encoding.

        Parameters
        ----------
        y : array-like of shape (n_samples,)
            Target values.

        Returns
        -------
        y_original : ndarray of shape (n_samples,)
            Original encoding.
        """
        check_is_fitted(self)
        xp, _ = get_namespace(y)
        y = column_or_1d(y, warn=True)
        # inverse transform of empty array is empty array
        if _num_samples(y) == 0:
            return xp.asarray([])

        diff = xpx.setdiff1d(
            y,
            xp.arange(self.classes_.shape[0], device=device(y)),
            xp=xp,
        )
        if diff.shape[0]:
            raise ValueError("y contains previously unseen labels: %s" % str(diff))
        y = xp.asarray(y)
        return xp.take(self.classes_, y, axis=0)

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.array_api_support = True
        tags.input_tags.two_d_array = False
        tags.target_tags.one_d_labels = True
        return tags


class LabelBinarizer(TransformerMixin, BaseEstimator, auto_wrap_output_keys=None):
    """Binarize labels in a one-vs-all fashion.

    Several regression and binary classification algorithms are
    available in scikit-learn. A simple way to extend these algorithms
    to the multi-class classification case is to use the so-called
    one-vs-all scheme.

    At learning time, this simply consists in learning one regressor
    or binary classifier per class. In doing so, one needs to convert
    multi-class labels to binary labels (belong or does not belong
    to the class). `LabelBinarizer` makes this process easy with the
    transform method.

    At prediction time, one assigns the class for which the corresponding
    model gave the greatest confidence. `LabelBinarizer` makes this easy
    with the :meth:`inverse_transform` method.

    Read more in the :ref:`User Guide <preprocessing_targets>`.

    Parameters
    ----------
    neg_label : int, default=0
        Value with which negative labels must be encoded.

    pos_label : int, default=1
        Value with which positive labels must be encoded.

    sparse_output : bool, default=False
        True if the returned array from transform is desired to be in sparse
        CSR format.

    Attributes
    ----------
    classes_ : ndarray of shape (n_classes,)
        Holds the label for each class.

    y_type_ : str
        Represents the type of the target data as evaluated by
        :func:`~sklearn.utils.multiclass.type_of_target`. Possible type are
        'continuous', 'continuous-multioutput', 'binary', 'multiclass',
        'multiclass-multioutput', 'multilabel-indicator', and 'unknown'.

    sparse_input_ : bool
        `True` if the input data to transform is given as a sparse matrix,
         `False` otherwise.

    See Also
    --------
    label_binarize : Function to perform the transform operation of
        LabelBinarizer with fixed classes.
    OneHotEncoder : Encode categorical features using a one-hot aka one-of-K
        scheme.

    Examples
    --------
    >>> from sklearn.preprocessing import LabelBinarizer
    >>> lb = LabelBinarizer()
    >>> lb.fit([1, 2, 6, 4, 2])
    LabelBinarizer()
    >>> lb.classes_
    array([1, 2, 4, 6])
    >>> lb.transform([1, 6])
    array([[1, 0, 0, 0],
           [0, 0, 0, 1]])

    Binary targets transform to a column vector

    >>> lb = LabelBinarizer()
    >>> lb.fit_transform(['yes', 'no', 'no', 'yes'])
    array([[1],
           [0],
           [0],
           [1]])

    Passing a 2D matrix for multilabel classification

    >>> import numpy as np
    >>> lb.fit(np.array([[0, 1, 1], [1, 0, 0]]))
    LabelBinarizer()
    >>> lb.classes_
    array([0, 1, 2])
    >>> lb.transform([0, 1, 2, 1])
    array([[1, 0, 0],
           [0, 1, 0],
           [0, 0, 1],
           [0, 1, 0]])
    """

    _parameter_constraints: dict = {
        "neg_label": [Integral],
        "pos_label": [Integral],
        "sparse_output": ["boolean"],
    }

    def __init__(self, *, neg_label=0, pos_label=1, sparse_output=False):
        self.neg_label = neg_label
        self.pos_label = pos_label
        self.sparse_output = sparse_output

    @_fit_context(prefer_skip_nested_validation=True)
    def fit(self, y):
        """Fit label binarizer.

        Parameters
        ----------
        y : ndarray of shape (n_samples,) or (n_samples, n_classes)
            Target values. The 2-d matrix should only contain 0 and 1,
            represents multilabel classification.

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        if self.neg_label >= self.pos_label:
            raise ValueError(
                f"neg_label={self.neg_label} must be strictly less than "
                f"pos_label={self.pos_label}."
            )

        if self.sparse_output and (self.pos_label == 0 or self.neg_label != 0):
            raise ValueError(
                "Sparse binarization is only supported with non "
                "zero pos_label and zero neg_label, got "
                f"pos_label={self.pos_label} and neg_label={self.neg_label}"
            )

        xp, is_array_api = get_namespace(y)

        if is_array_api and self.sparse_output and not _is_numpy_namespace(xp):
            raise ValueError(
                "`sparse_output=True` is not supported for array API "
                f"namespace {xp.__name__}. "
                "Use `sparse_output=False` to return a dense array instead."
            )

        self.y_type_ = type_of_target(y, input_name="y")

        if "multioutput" in self.y_type_:
            raise ValueError(
                "Multioutput target data is not supported with label binarization"
            )
        if _num_samples(y) == 0:
            raise ValueError("y has 0 samples: %r" % y)

        self.sparse_input_ = sp.issparse(y)
        self.classes_ = unique_labels(y)
        return self

    def fit_transform(self, y):
        """Fit label binarizer/transform multi-class labels to binary labels.

        The output of transform is sometimes referred to as
        the 1-of-K coding scheme.

        Parameters
        ----------
        y : {ndarray, sparse matrix} of shape (n_samples,) or \
                (n_samples, n_classes)
            Target values. The 2-d matrix should only contain 0 and 1,
            represents multilabel classification. Sparse matrix can be
            CSR, CSC, COO, DOK, or LIL.

        Returns
        -------
        Y : {ndarray, sparse matrix} of shape (n_samples, n_classes)
            Shape will be (n_samples, 1) for binary problems. Sparse matrix
            will be of CSR format.
        """
        return self.fit(y).transform(y)

    def transform(self, y):
        """Transform multi-class labels to binary labels.

        The output of transform is sometimes referred to by some authors as
        the 1-of-K coding scheme.

        Parameters
        ----------
        y : {array, sparse matrix} of shape (n_samples,) or \
                (n_samples, n_classes)
            Target values. The 2-d matrix should only contain 0 and 1,
            represents multilabel classification. Sparse matrix can be
            CSR, CSC, COO, DOK, or LIL.

        Returns
        -------
        Y : {ndarray, sparse matrix} of shape (n_samples, n_classes)
            Shape will be (n_samples, 1) for binary problems. Sparse matrix
            will be of CSR format.
        """
        check_is_fitted(self)

        xp, is_array_api = get_namespace(y)

        if is_array_api and self.sparse_output and not _is_numpy_namespace(xp):
            raise ValueError(
                "`sparse_output=True` is not supported for array API "
                f"namespace {xp.__name__}. "
                "Use `sparse_output=False` to return a dense array instead."
            )

        y_is_multilabel = type_of_target(y).startswith("multilabel")
        if y_is_multilabel and not self.y_type_.startswith("multilabel"):
            raise ValueError("The object was not fitted with multilabel input.")

        return label_binarize(
            y,
            classes=self.classes_,
            pos_label=self.pos_label,
            neg_label=self.neg_label,
            sparse_output=self.sparse_output,
        )

    def inverse_transform(self, Y, threshold=None):
        """Transform binary labels back to multi-class labels.

        Parameters
        ----------
        Y : {ndarray, sparse matrix} of shape (n_samples, n_classes)
            Target values. All sparse matrices are converted to CSR before
            inverse transformation.

        threshold : float, default=None
            Threshold used in the binary and multi-label cases.

            Use 0 when ``Y`` contains the output of :term:`decision_function`
            (classifier).
            Use 0.5 when ``Y`` contains the output of :term:`predict_proba`.

            If None, the threshold is assumed to be half way between
            neg_label and pos_label.

        Returns
        -------
        y_original : {ndarray, sparse matrix} of shape (n_samples,)
            Target values. Sparse matrix will be of CSR format.

        Notes
        -----
        In the case when the binary labels are fractional
        (probabilistic), :meth:`inverse_transform` chooses the class with the
        greatest value. Typically, this allows to use the output of a
        linear model's :term:`decision_function` method directly as the input
        of :meth:`inverse_transform`.
        """
        check_is_fitted(self)

        xp, is_array_api = get_namespace(Y)

        if is_array_api and self.sparse_input_ and not _is_numpy_namespace(xp):
            raise ValueError(
                "`LabelBinarizer` was fitted on a sparse matrix, and therefore cannot "
                f"inverse transform a {xp.__name__} array back to a sparse matrix."
            )

        if threshold is None:
            threshold = (self.pos_label + self.neg_label) / 2.0

        if self.y_type_ == "multiclass":
            y_inv = _inverse_binarize_multiclass(Y, self.classes_, xp=xp)
        else:
            y_inv = _inverse_binarize_thresholding(
                Y, self.y_type_, self.classes_, threshold, xp=xp
            )

        if self.sparse_input_:
            y_inv = sp.csr_matrix(y_inv)
        elif sp.issparse(y_inv):
            y_inv = y_inv.toarray()

        return y_inv

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.input_tags.two_d_array = False
        tags.target_tags.one_d_labels = True
        return tags


@validate_params(
    {
        "y": ["array-like", "sparse matrix"],
        "classes": ["array-like"],
        "neg_label": [Interval(Integral, None, None, closed="neither")],
        "pos_label": [Interval(Integral, None, None, closed="neither")],
        "sparse_output": ["boolean"],
    },
    prefer_skip_nested_validation=True,
)
def label_binarize(y, *, classes, neg_label=0, pos_label=1, sparse_output=False):
    """Binarize labels in a one-vs-all fashion.

    Several regression and binary classification algorithms are
    available in scikit-learn. A simple way to extend these algorithms
    to the multi-class classification case is to use the so-called
    one-vs-all scheme.

    This function makes it possible to compute this transformation for a
    fixed set of class labels known ahead of time.

    Parameters
    ----------
    y : array-like or sparse matrix
        Sequence of integer labels or multilabel data to encode.

    classes : array-like of shape (n_classes,)
        Uniquely holds the label for each class.

    neg_label : int, default=0
        Value with which negative labels must be encoded.

    pos_label : int, default=1
        Value with which positive labels must be encoded.

    sparse_output : bool, default=False,
        Set to true if output binary array is desired in CSR sparse format.

    Returns
    -------
    Y : {ndarray, sparse matrix} of shape (n_samples, n_classes)
        Shape will be (n_samples, 1) for binary problems. Sparse matrix will
        be of CSR format.

    See Also
    --------
    LabelBinarizer : Class used to wrap the functionality of label_binarize and
        allow for fitting to classes independently of the transform operation.

    Examples
    --------
    >>> from sklearn.preprocessing import label_binarize
    >>> label_binarize([1, 6], classes=[1, 2, 4, 6])
    array([[1, 0, 0, 0],
           [0, 0, 0, 1]])

    The class ordering is preserved:

    >>> label_binarize([1, 6], classes=[1, 6, 4, 2])
    array([[1, 0, 0, 0],
           [0, 1, 0, 0]])

    Binary targets transform to a column vector

    >>> label_binarize(['yes', 'no', 'no', 'yes'], classes=['no', 'yes'])
    array([[1],
           [0],
           [0],
           [1]])
    """
    if not isinstance(y, list):
        # XXX Workaround that will be removed when list of list format is
        # dropped
        y = check_array(
            y, input_name="y", accept_sparse="csr", ensure_2d=False, dtype=None
        )
    else:
        if _num_samples(y) == 0:
            raise ValueError("y has 0 samples: %r" % y)
    if neg_label >= pos_label:
        raise ValueError(
            "neg_label={0} must be strictly less than pos_label={1}.".format(
                neg_label, pos_label
            )
        )

    if sparse_output and (pos_label == 0 or neg_label != 0):
        raise ValueError(
            "Sparse binarization is only supported with non "
            "zero pos_label and zero neg_label, got "
            "pos_label={0} and neg_label={1}"
            "".format(pos_label, neg_label)
        )

    # To account for pos_label == 0 in the dense case
    pos_switch = pos_label == 0
    if pos_switch:
        pos_label = -neg_label

    y_type = type_of_target(y)
    if "multioutput" in y_type:
        raise ValueError(
            "Multioutput target data is not supported with label binarization"
        )
    if y_type == "unknown":
        raise ValueError("The type of target data is not known")

    xp, is_array_api, device_ = get_namespace_and_device(y)

    if is_array_api and sparse_output and not _is_numpy_namespace(xp):
        raise ValueError(
            "`sparse_output=True` is not supported for array API "
            f"'namespace {xp.__name__}'. "
            "Use `sparse_output=False` to return a dense array instead."
        )

    try:
        classes = xp.asarray(classes, device=device_)
    except (ValueError, TypeError) as e:
        # `classes` contains an unsupported dtype for this namespace.
        # For example, attempting to create torch.tensor(["yes", "no"]) will fail.
        raise ValueError(
            f"`classes` contains unsupported dtype for array API namespace "
            f"'{xp.__name__}'."
        ) from e

    n_samples = y.shape[0] if hasattr(y, "shape") else len(y)
    n_classes = classes.shape[0]
    if hasattr(y, "dtype") and xp.isdtype(y.dtype, "integral"):
        int_dtype_ = y.dtype
    else:
        int_dtype_ = indexing_dtype(xp)

    if y_type == "binary":
        if n_classes == 1:
            if sparse_output:
                return sp.csr_matrix((n_samples, 1), dtype=int)
            else:
                Y = xp.zeros((n_samples, 1), dtype=int_dtype_)
                Y += neg_label
                return Y
        elif n_classes >= 3:
            y_type = "multiclass"

    sorted_class = xp.sort(classes)
    if y_type == "multilabel-indicator":
        y_n_classes = y.shape[1] if hasattr(y, "shape") else len(y[0])
        if n_classes != y_n_classes:
            raise ValueError(
                "classes {0} mismatch with the labels {1} found in the data".format(
                    classes, unique_labels(y)
                )
            )

    if y_type in ("binary", "multiclass"):
        y = column_or_1d(y)

        # pick out the known labels from y
        y_in_classes = _isin(y, classes, xp=xp)
        y_seen = y[y_in_classes]
        indices = xp.searchsorted(sorted_class, y_seen)
        # cast `y_in_classes` to integer dtype for `xp.cumulative_sum`
        y_in_classes = xp.astype(y_in_classes, int_dtype_)
        indptr = xp.concat(
            (
                xp.asarray([0], device=device_),
                xp.cumulative_sum(y_in_classes, axis=0),
            )
        )
        data = xp.full_like(indices, pos_label)

        # Use NumPy to construct the sparse matrix of one-hot labels
        Y = sp.csr_matrix(
            (
                _convert_to_numpy(data, xp=xp),
                _convert_to_numpy(indices, xp=xp),
                _convert_to_numpy(indptr, xp=xp),
            ),
            shape=(n_samples, n_classes),
        )

        if not sparse_output:
            Y = xp.asarray(Y.toarray(), device=device_)

    elif y_type == "multilabel-indicator":
        if sparse_output:
            Y = sp.csr_matrix(y)
            if pos_label != 1:
                data = xp.full_like(Y.data, pos_label)
                Y.data = data
        else:
            if sp.issparse(y):
                y = y.toarray()

            Y = xp.asarray(y, device=device_, copy=True)
            if pos_label != 1:
                Y[Y != 0] = pos_label

    else:
        raise ValueError(
            "%s target data is not supported with label binarization" % y_type
        )

    if not sparse_output:
        if neg_label != 0:
            Y[Y == 0] = neg_label

        if pos_switch:
            Y[Y == pos_label] = 0

        Y = xp.astype(Y, int_dtype_, copy=False)
    else:
        Y.data = Y.data.astype(int, copy=False)

    # preserve label ordering
    if xp.any(classes != sorted_class):
        indices = xp.searchsorted(sorted_class, classes)
        Y = Y[:, indices]

    if y_type == "binary":
        if sparse_output:
            Y = Y[:, [-1]]
        else:
            Y = xp.reshape(Y[:, -1], (-1, 1))

    return Y


def _inverse_binarize_multiclass(y, classes, xp=None):
    """Inverse label binarization transformation for multiclass.

    Multiclass uses the maximal score instead of a threshold.
    """
    if sp.issparse(y):
        classes = np.asarray(classes)
        # Find the argmax for each row in y where y is a CSR matrix

        y = y.tocsr()
        n_samples, n_outputs = y.shape
        outputs = np.arange(n_outputs)
        row_max = min_max_axis(y, 1)[1]
        row_nnz = np.diff(y.indptr)

        y_data_repeated_max = np.repeat(row_max, row_nnz)
        # picks out all indices obtaining the maximum per row
        y_i_all_argmax = np.flatnonzero(y_data_repeated_max == y.data)

        # For corner case where last row has a max of 0
        if row_max[-1] == 0:
            y_i_all_argmax = np.append(y_i_all_argmax, [len(y.data)])

        # Gets the index of the first argmax in each row from y_i_all_argmax
        index_first_argmax = np.searchsorted(y_i_all_argmax, y.indptr[:-1])
        # first argmax of each row
        y_ind_ext = np.append(y.indices, [0])
        y_i_argmax = y_ind_ext[y_i_all_argmax[index_first_argmax]]
        # Handle rows of all 0
        y_i_argmax[np.where(row_nnz == 0)[0]] = 0

        # Handles rows with max of 0 that contain negative numbers
        samples = np.arange(n_samples)[(row_nnz > 0) & (row_max.ravel() == 0)]
        for i in samples:
            ind = y.indices[y.indptr[i] : y.indptr[i + 1]]
            y_i_argmax[i] = classes[np.setdiff1d(outputs, ind)][0]

        return classes[y_i_argmax]
    else:
        xp, _, device_ = get_namespace_and_device(y, xp=xp)
        classes = xp.asarray(classes, device=device_)
        indices = xp.argmax(y, axis=1)
        indices = xp.clip(indices, 0, classes.shape[0] - 1)

        return classes[indices]


def _inverse_binarize_thresholding(y, output_type, classes, threshold, xp=None):
    """Inverse label binarization transformation using thresholding."""

    if output_type == "binary" and y.ndim == 2 and y.shape[1] > 2:
        raise ValueError("output_type='binary', but y.shape = {0}".format(y.shape))

    xp, _, device_ = get_namespace_and_device(y, xp=xp)
    classes = xp.asarray(classes, device=device_)

    if output_type != "binary" and y.shape[1] != classes.shape[0]:
        raise ValueError(
            "The number of class is not equal to the number of dimension of y."
        )

    dtype_ = _find_matching_floating_dtype(y, xp=xp)
    if hasattr(y, "dtype") and xp.isdtype(y.dtype, "integral"):
        int_dtype_ = y.dtype
    else:
        int_dtype_ = indexing_dtype(xp)

    # Perform thresholding
    if sp.issparse(y):
        if threshold > 0:
            if y.format not in ("csr", "csc"):
                y = y.tocsr()
            y.data = np.array(y.data > threshold, dtype=int)
            y.eliminate_zeros()
        else:
            y = xp.asarray(y.toarray() > threshold, dtype=int_dtype_, device=device_)
    else:
        y = xp.asarray(
            xp.asarray(y, dtype=dtype_, device=device_) > threshold,
            dtype=int_dtype_,
            device=device_,
        )

    # Inverse transform data
    if output_type == "binary":
        if sp.issparse(y):
            y = y.toarray()
        if y.ndim == 2 and y.shape[1] == 2:
            return classes[y[:, 1]]
        else:
            if classes.shape[0] == 1:
                return xp.repeat(classes[0], len(y))
            else:
                return classes[xp.reshape(y, (-1,))]

    elif output_type == "multilabel-indicator":
        return y

    else:
        raise ValueError("{0} format is not supported".format(output_type))


class MultiLabelBinarizer(TransformerMixin, BaseEstimator, auto_wrap_output_keys=None):
    """Transform between iterable of iterables and a multilabel format.

    Although a list of sets or tuples is a very intuitive format for multilabel
    data, it is unwieldy to process. This transformer converts between this
    intuitive format and the supported multilabel format: a (samples x classes)
    binary matrix indicating the presence of a class label.

    Parameters
    ----------
    classes : array-like of shape (n_classes,), default=None
        Indicates an ordering for the class labels.
        All entries should be unique (cannot contain duplicate classes).

    sparse_output : bool, default=False
        Set to True if output binary array is desired in CSR sparse format.

    Attributes
    ----------
    classes_ : ndarray of shape (n_classes,)
        A copy of the `classes` parameter when provided.
        Otherwise it corresponds to the sorted set of classes found
        when fitting.

    See Also
    --------
    OneHotEncoder : Encode categorical features using a one-hot aka one-of-K
        scheme.

    Examples
    --------
    >>> from sklearn.preprocessing import MultiLabelBinarizer
    >>> mlb = MultiLabelBinarizer()
    >>> mlb.fit_transform([(1, 2), (3,)])
    array([[1, 1, 0],
           [0, 0, 1]])
    >>> mlb.classes_
    array([1, 2, 3])

    >>> mlb.fit_transform([{'sci-fi', 'thriller'}, {'comedy'}])
    array([[0, 1, 1],
           [1, 0, 0]])
    >>> list(mlb.classes_)
    ['comedy', 'sci-fi', 'thriller']

    A common mistake is to pass in a list, which leads to the following issue:

    >>> mlb = MultiLabelBinarizer()
    >>> mlb.fit(['sci-fi', 'thriller', 'comedy'])
    MultiLabelBinarizer()
    >>> mlb.classes_
    array(['-', 'c', 'd', 'e', 'f', 'h', 'i', 'l', 'm', 'o', 'r', 's', 't',
        'y'], dtype=object)

    To correct this, the list of labels should be passed in as:

    >>> mlb = MultiLabelBinarizer()
    >>> mlb.fit([['sci-fi', 'thriller', 'comedy']])
    MultiLabelBinarizer()
    >>> mlb.classes_
    array(['comedy', 'sci-fi', 'thriller'], dtype=object)
    """

    _parameter_constraints: dict = {
        "classes": ["array-like", None],
        "sparse_output": ["boolean"],
    }

    def __init__(self, *, classes=None, sparse_output=False):
        self.classes = classes
        self.sparse_output = sparse_output

    @_fit_context(prefer_skip_nested_validation=True)
    def fit(self, y):
        """Fit the label sets binarizer, storing :term:`classes_`.

        Parameters
        ----------
        y : iterable of iterables
            A set of labels (any orderable and hashable object) for each
            sample. If the `classes` parameter is set, `y` will not be
            iterated.

        Returns
        -------
        self : object
            Fitted estimator.
        """
        self._cached_dict = None

        if self.classes is None:
            classes = sorted(set(itertools.chain.from_iterable(y)))
        elif len(set(self.classes)) < len(self.classes):
            raise ValueError(
                "The classes argument contains duplicate "
                "classes. Remove these duplicates before passing "
                "them to MultiLabelBinarizer."
            )
        else:
            classes = self.classes
        dtype = int if all(isinstance(c, int) for c in classes) else object
        self.classes_ = np.empty(len(classes), dtype=dtype)
        self.classes_[:] = classes
        return self

    @_fit_context(prefer_skip_nested_validation=True)
    def fit_transform(self, y):
        """Fit the label sets binarizer and transform the given label sets.

        Parameters
        ----------
        y : iterable of iterables
            A set of labels (any orderable and hashable object) for each
            sample. If the `classes` parameter is set, `y` will not be
            iterated.

        Returns
        -------
        y_indicator : {ndarray, sparse matrix} of shape (n_samples, n_classes)
            A matrix such that `y_indicator[i, j] = 1` iff `classes_[j]`
            is in `y[i]`, and 0 otherwise. Sparse matrix will be of CSR
            format.
        """
        if self.classes is not None:
            return self.fit(y).transform(y)

        self._cached_dict = None

        # Automatically increment on new class
        class_mapping = defaultdict(int)
        class_mapping.default_factory = class_mapping.__len__
        yt = self._transform(y, class_mapping)

        # sort classes and reorder columns
        tmp = sorted(class_mapping, key=class_mapping.get)

        # (make safe for tuples)
        dtype = int if all(isinstance(c, int) for c in tmp) else object
        class_mapping = np.empty(len(tmp), dtype=dtype)
        class_mapping[:] = tmp
        self.classes_, inverse = np.unique(class_mapping, return_inverse=True)
        # ensure yt.indices keeps its current dtype
        yt.indices = np.asarray(inverse[yt.indices], dtype=yt.indices.dtype)

        if not self.sparse_output:
            yt = yt.toarray()

        return yt

    def transform(self, y):
        """Transform the given label sets.

        Parameters
        ----------
        y : iterable of iterables
            A set of labels (any orderable and hashable object) for each
            sample. If the `classes` parameter is set, `y` will not be
            iterated.

        Returns
        -------
        y_indicator : array or CSR matrix, shape (n_samples, n_classes)
            A matrix such that `y_indicator[i, j] = 1` iff `classes_[j]` is in
            `y[i]`, and 0 otherwise.
        """
        check_is_fitted(self)

        class_to_index = self._build_cache()
        yt = self._transform(y, class_to_index)

        if not self.sparse_output:
            yt = yt.toarray()

        return yt

    def _build_cache(self):
        if self._cached_dict is None:
            self._cached_dict = dict(zip(self.classes_, range(len(self.classes_))))

        return self._cached_dict

    def _transform(self, y, class_mapping):
        """Transforms the label sets with a given mapping.

        Parameters
        ----------
        y : iterable of iterables
            A set of labels (any orderable and hashable object) for each
            sample. If the `classes` parameter is set, `y` will not be
            iterated.

        class_mapping : Mapping
            Maps from label to column index in label indicator matrix.

        Returns
        -------
        y_indicator : sparse matrix of shape (n_samples, n_classes)
            Label indicator matrix. Will be of CSR format.
        """
        indices = array.array("i")
        indptr = array.array("i", [0])
        unknown = set()
        for labels in y:
            index = set()
            for label in labels:
                try:
                    index.add(class_mapping[label])
                except KeyError:
                    unknown.add(label)
            indices.extend(index)
            indptr.append(len(indices))
        if unknown:
            warnings.warn(
                "unknown class(es) {0} will be ignored".format(sorted(unknown, key=str))
            )
        data = np.ones(len(indices), dtype=int)

        return sp.csr_matrix(
            (data, indices, indptr), shape=(len(indptr) - 1, len(class_mapping))
        )

    def inverse_transform(self, yt):
        """Transform the given indicator matrix into label sets.

        Parameters
        ----------
        yt : {ndarray, sparse matrix} of shape (n_samples, n_classes)
            A matrix containing only 1s ands 0s.

        Returns
        -------
        y_original : list of tuples
            The set of labels for each sample such that `y[i]` consists of
            `classes_[j]` for each `yt[i, j] == 1`.
        """
        check_is_fitted(self)

        if yt.shape[1] != len(self.classes_):
            raise ValueError(
                "Expected indicator for {0} classes, but got {1}".format(
                    len(self.classes_), yt.shape[1]
                )
            )

        if sp.issparse(yt):
            yt = yt.tocsr()
            if len(yt.data) != 0 and len(np.setdiff1d(yt.data, [0, 1])) > 0:
                raise ValueError("Expected only 0s and 1s in label indicator.")
            return [
                tuple(self.classes_.take(yt.indices[start:end]))
                for start, end in zip(yt.indptr[:-1], yt.indptr[1:])
            ]
        else:
            unexpected = np.setdiff1d(yt, [0, 1])
            if len(unexpected) > 0:
                raise ValueError(
                    "Expected only 0s and 1s in label indicator. Also got {0}".format(
                        unexpected
                    )
                )
            return [tuple(self.classes_.compress(indicators)) for indicators in yt]

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.input_tags.two_d_array = False
        tags.target_tags.two_d_labels = True
        return tags
