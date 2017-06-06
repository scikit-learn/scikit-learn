# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Mathieu Blondel <mathieu@mblondel.org>
#          Olivier Grisel <olivier.grisel@ensta.org>
#          Andreas Mueller <amueller@ais.uni-bonn.de>
#          Joel Nothman <joel.nothman@gmail.com>
#          Hamzeh Alsalhi <ha258@cornell.edu>
# License: BSD 3 clause

from collections import defaultdict
import itertools
import array

import numpy as np
import scipy.sparse as sp

from ..base import (BaseEstimator, TransformerMixin, RegressorMixin, clone,
                    is_regressor)

from ..utils.fixes import np_version
from ..utils.fixes import sparse_min_max
from ..utils.fixes import astype
from ..utils.fixes import in1d
from ..utils import column_or_1d
from ..utils.validation import check_array
from ..utils.validation import check_is_fitted
from ..utils.validation import _num_samples
from ..utils.multiclass import unique_labels
from ..utils.multiclass import type_of_target

from ..externals import six

from ._function_transformer import FunctionTransformer

zip = six.moves.zip
map = six.moves.map

__all__ = [
    'label_binarize',
    'LabelBinarizer',
    'LabelEncoder',
    'MultiLabelBinarizer',
    'TransformedTargetRegressor'
]


def _check_numpy_unicode_bug(labels):
    """Check that user is not subject to an old numpy bug

    Fixed in master before 1.7.0:

      https://github.com/numpy/numpy/pull/243

    """
    if np_version[:3] < (1, 7, 0) and labels.dtype.kind == 'U':
        raise RuntimeError("NumPy < 1.7.0 does not implement searchsorted"
                           " on unicode data correctly. Please upgrade"
                           " NumPy to use LabelEncoder with unicode inputs.")


class TransformedTargetRegressor(BaseEstimator, RegressorMixin):
    """Meta-estimator to apply a transformation to the target before fitting.

    Useful for applying a non-linear transformation in regression
    problems. This transformation can be given as a Transformer such as the
    QuantileTransformer or as a function and its inverse such as ``np.log`` and
    ``np.exp``.

    The computation during ``fit`` is::

        estimator.fit(X, func(y))

    or::

        estimator.fit(X, transformer.transform(y))

    The computation during ``predict`` is::

        inverse_func(estimator.predict(X))

    or::

        transformer.inverse_transform(estimator.predict(X))

    Parameters
    ----------
    estimator : object, (default=LinearRegression())
        Estimator object derived from ``RegressorMixin``.

    transformer : object, (default=None)
        Estimator object derived from ``TransformerMixin``. Cannot be set at
        the same time as ``func`` and ``inverse_func``. If ``None`` and
        ``func`` and ``inverse_func`` as well, the transformer will be an
        identity transformer.

    func : function, (default=None)
        Function to apply to ``y`` before passing to ``fit``. Cannot be set at
        the same time than ``transformer``. If ``None`` and ``transformer`` as
        well, the function used will be the identity function.

    inverse_func : function, (default=None)
        Function apply to the prediction of the estimator. Cannot be set at
        the same time than ``transformer``. If ``None`` and ``transformer`` as
        well, the function used will be the identity function.

    check_invertible : bool, (default=True)
        Whether to check that ``transform`` followed by ``inverse_transform``
        or ``func`` followed by ``inverse_func`` lead to the original data.

    Attributes
    ----------
    estimator_ : object
        Fitted estimator.

    transformer_ : object
        Used transformer in ``fit`` and ``predict``.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.linear_model import LinearRegression
    >>> from sklearn.preprocessing.label import TransformedTargetRegressor
    >>> tt = TransformedTargetRegressor(estimator=LinearRegression(),
    ...                                 func=np.log, inverse_func=np.exp)
    >>> X = np.arange(4).reshape(-1, 1)
    >>> y = np.exp(2 * X).ravel()
    >>> tt.fit(X, y)
    ... #doctest: +NORMALIZE_WHITESPACE
    TransformedTargetRegressor(check_invertible=True,
        estimator=LinearRegression(copy_X=True, fit_intercept=True, n_jobs=1,
                                   normalize=False),
        func=<ufunc 'log'>, inverse_func=<ufunc 'exp'>, transformer=None)
    >>> tt.score(X, y)
    1.0
    >>> tt.estimator_.coef_
    array([ 2.])

    """
    def __init__(self, estimator, transformer=None,
                 func=None, inverse_func=None, check_invertible=True):
        self.estimator = estimator
        self.transformer = transformer
        self.func = func
        self.inverse_func = inverse_func
        self.check_invertible = check_invertible
        # we probably need to change this ones we have tags
        self._estimator_type = estimator._estimator_type

    def _validate_transformer(self, y):
        if (self.transformer is not None and
                (self.func is not None or self.inverse_func is not None)):
            raise ValueError("Both 'transformer' and functions 'func'/"
                             "'inverse_func' cannot be set at the same time.")
        elif self.transformer is not None:
            self.transformer_ = clone(self.transformer)
        else:
            self.transformer_ = FunctionTransformer(
                func=self.func, inverse_func=self.inverse_func, validate=False)
        self.transformer_.fit(y)
        # XXX: is it a necessary test? one might not want an invertible
        # function.
        if self.check_invertible:
            n_subsample = min(1000, y.shape[0])
            subsample_idx = np.random.randint(y.shape[0], size=n_subsample)
            diff = np.abs((y[subsample_idx] -
                           self.transformer_.inverse_transform(
                               self.transformer_.transform(y[subsample_idx]))))
            if np.sum(diff) > 1e-7:
                raise ValueError("The provided functions or transformer are"
                                 " not strictly invertible.")

    def fit(self, X, y, sample_weight=None):
        """Fit the model according to the given training data.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training vector, where n_samples is the number of samples and
            n_features is the number of features.

        y : array-like, shape (n_samples,)
            Target vector relative to X.

        sample_weight : array-like, shape (n_samples,) optional
            Array of weights that are assigned to individual samples.
            If not provided, then each sample is given unit weight.

        Returns
        -------
        self : object
            Returns self.
        """
        self._validate_transformer(y)
        self.estimator_ = clone(self.estimator)
        self.estimator_.fit(X, self.transformer_.transform(y),
                            sample_weight=sample_weight)
        return self

    def predict(self, X):
        """Predict using the base estimator, applying inverse.

        The estimator is used to predict and the ``inverse_func`` or
        ``inverse_transform`` is applied before returning the prediction.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = (n_samples, n_features)
            Samples.

        Returns
        -------
        y_hat : array, shape = (n_samples,)
            Predicted values.

        """
        check_is_fitted(self, "estimator_")
        return self.transformer_.inverse_transform(self.estimator_.predict(X))

    def score(self, X, y, sample_weight=None):
        """Returns the coefficient of determination R^2 of the prediction.

        The coefficient R^2 is defined as (1 - u/v), where u is the regression
        sum of squares ((y_true - y_pred) ** 2).sum() and v is the residual sum
        of squares ((y_true - y_true.mean()) ** 2).sum().  Best possible score
        is 1.0 and it can be negative (because the model can be arbitrarily
        worse). A constant model that always predicts the expected value of y,
        disregarding the input features, would get a R^2 score of 0.0.

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            Test samples.

        y : array-like, shape = (n_samples) or (n_samples, n_outputs)
            True values for X.

        sample_weight : array-like, shape = [n_samples], optional
            Sample weights.

        Returns
        -------
        score : float
            R^2 of self.predict(X) wrt. y.

        """

        check_is_fitted(self, "estimator_")
        if not is_regressor(self.estimator_):
            if not hasattr(self.estimator_, "_estimator_type"):
                err = "estimator has declared no _estimator_type."
            else:
                err = "estimator has _estimator_type {}".format(
                    self.estimator_._estimator_type)
            raise NotImplementedError("TransformedTargetRegressor should be a"
                                      " regressor. This " + err)
        else:
            return super(TransformedTargetRegressor, self).score(X, y,
                                                                 sample_weight)


class LabelEncoder(BaseEstimator, TransformerMixin):
    """Encode labels with value between 0 and n_classes-1.

    Read more in the :ref:`User Guide <preprocessing_targets>`.

    Attributes
    ----------
    classes_ : array of shape (n_class,)
        Holds the label for each class.

    Examples
    --------
    `LabelEncoder` can be used to normalize labels.

    >>> from sklearn import preprocessing
    >>> le = preprocessing.LabelEncoder()
    >>> le.fit([1, 2, 2, 6])
    LabelEncoder()
    >>> le.classes_
    array([1, 2, 6])
    >>> le.transform([1, 1, 2, 6]) #doctest: +ELLIPSIS
    array([0, 0, 1, 2]...)
    >>> le.inverse_transform([0, 0, 1, 2])
    array([1, 1, 2, 6])

    It can also be used to transform non-numerical labels (as long as they are
    hashable and comparable) to numerical labels.

    >>> le = preprocessing.LabelEncoder()
    >>> le.fit(["paris", "paris", "tokyo", "amsterdam"])
    LabelEncoder()
    >>> list(le.classes_)
    ['amsterdam', 'paris', 'tokyo']
    >>> le.transform(["tokyo", "tokyo", "paris"]) #doctest: +ELLIPSIS
    array([2, 2, 1]...)
    >>> list(le.inverse_transform([2, 2, 1]))
    ['tokyo', 'tokyo', 'paris']

    See also
    --------
    sklearn.preprocessing.OneHotEncoder : encode categorical integer features
        using a one-hot aka one-of-K scheme.
    """

    def fit(self, y):
        """Fit label encoder

        Parameters
        ----------
        y : array-like of shape (n_samples,)
            Target values.

        Returns
        -------
        self : returns an instance of self.
        """
        y = column_or_1d(y, warn=True)
        _check_numpy_unicode_bug(y)
        self.classes_ = np.unique(y)
        return self

    def fit_transform(self, y):
        """Fit label encoder and return encoded labels

        Parameters
        ----------
        y : array-like of shape [n_samples]
            Target values.

        Returns
        -------
        y : array-like of shape [n_samples]
        """
        y = column_or_1d(y, warn=True)
        _check_numpy_unicode_bug(y)
        self.classes_, y = np.unique(y, return_inverse=True)
        return y

    def transform(self, y):
        """Transform labels to normalized encoding.

        Parameters
        ----------
        y : array-like of shape [n_samples]
            Target values.

        Returns
        -------
        y : array-like of shape [n_samples]
        """
        check_is_fitted(self, 'classes_')
        y = column_or_1d(y, warn=True)

        classes = np.unique(y)
        _check_numpy_unicode_bug(classes)
        if len(np.intersect1d(classes, self.classes_)) < len(classes):
            diff = np.setdiff1d(classes, self.classes_)
            raise ValueError("y contains new labels: %s" % str(diff))
        return np.searchsorted(self.classes_, y)

    def inverse_transform(self, y):
        """Transform labels back to original encoding.

        Parameters
        ----------
        y : numpy array of shape [n_samples]
            Target values.

        Returns
        -------
        y : numpy array of shape [n_samples]
        """
        check_is_fitted(self, 'classes_')

        diff = np.setdiff1d(y, np.arange(len(self.classes_)))
        if diff:
            raise ValueError("y contains new labels: %s" % str(diff))
        y = np.asarray(y)
        return self.classes_[y]


class LabelBinarizer(BaseEstimator, TransformerMixin):
    """Binarize labels in a one-vs-all fashion

    Several regression and binary classification algorithms are
    available in the scikit. A simple way to extend these algorithms
    to the multi-class classification case is to use the so-called
    one-vs-all scheme.

    At learning time, this simply consists in learning one regressor
    or binary classifier per class. In doing so, one needs to convert
    multi-class labels to binary labels (belong or does not belong
    to the class). LabelBinarizer makes this process easy with the
    transform method.

    At prediction time, one assigns the class for which the corresponding
    model gave the greatest confidence. LabelBinarizer makes this easy
    with the inverse_transform method.

    Read more in the :ref:`User Guide <preprocessing_targets>`.

    Parameters
    ----------

    neg_label : int (default: 0)
        Value with which negative labels must be encoded.

    pos_label : int (default: 1)
        Value with which positive labels must be encoded.

    sparse_output : boolean (default: False)
        True if the returned array from transform is desired to be in sparse
        CSR format.

    Attributes
    ----------

    classes_ : array of shape [n_class]
        Holds the label for each class.

    y_type_ : str,
        Represents the type of the target data as evaluated by
        utils.multiclass.type_of_target. Possible type are 'continuous',
        'continuous-multioutput', 'binary', 'multiclass',
        'multiclass-multioutput', 'multilabel-indicator', and 'unknown'.

    sparse_input_ : boolean,
        True if the input data to transform is given as a sparse matrix, False
        otherwise.

    Examples
    --------
    >>> from sklearn import preprocessing
    >>> lb = preprocessing.LabelBinarizer()
    >>> lb.fit([1, 2, 6, 4, 2])
    LabelBinarizer(neg_label=0, pos_label=1, sparse_output=False)
    >>> lb.classes_
    array([1, 2, 4, 6])
    >>> lb.transform([1, 6])
    array([[1, 0, 0, 0],
           [0, 0, 0, 1]])

    Binary targets transform to a column vector

    >>> lb = preprocessing.LabelBinarizer()
    >>> lb.fit_transform(['yes', 'no', 'no', 'yes'])
    array([[1],
           [0],
           [0],
           [1]])

    Passing a 2D matrix for multilabel classification

    >>> import numpy as np
    >>> lb.fit(np.array([[0, 1, 1], [1, 0, 0]]))
    LabelBinarizer(neg_label=0, pos_label=1, sparse_output=False)
    >>> lb.classes_
    array([0, 1, 2])
    >>> lb.transform([0, 1, 2, 1])
    array([[1, 0, 0],
           [0, 1, 0],
           [0, 0, 1],
           [0, 1, 0]])

    See also
    --------
    label_binarize : function to perform the transform operation of
        LabelBinarizer with fixed classes.
    sklearn.preprocessing.OneHotEncoder : encode categorical integer features
        using a one-hot aka one-of-K scheme.
    """

    def __init__(self, neg_label=0, pos_label=1, sparse_output=False):
        if neg_label >= pos_label:
            raise ValueError("neg_label={0} must be strictly less than "
                             "pos_label={1}.".format(neg_label, pos_label))

        if sparse_output and (pos_label == 0 or neg_label != 0):
            raise ValueError("Sparse binarization is only supported with non "
                             "zero pos_label and zero neg_label, got "
                             "pos_label={0} and neg_label={1}"
                             "".format(pos_label, neg_label))

        self.neg_label = neg_label
        self.pos_label = pos_label
        self.sparse_output = sparse_output

    def fit(self, y):
        """Fit label binarizer

        Parameters
        ----------
        y : array of shape [n_samples,] or [n_samples, n_classes]
            Target values. The 2-d matrix should only contain 0 and 1,
            represents multilabel classification.

        Returns
        -------
        self : returns an instance of self.
        """
        self.y_type_ = type_of_target(y)
        if 'multioutput' in self.y_type_:
            raise ValueError("Multioutput target data is not supported with "
                             "label binarization")
        if _num_samples(y) == 0:
            raise ValueError('y has 0 samples: %r' % y)

        self.sparse_input_ = sp.issparse(y)
        self.classes_ = unique_labels(y)
        return self

    def fit_transform(self, y):
        """Fit label binarizer and transform multi-class labels to binary
        labels.

        The output of transform is sometimes referred to    as
        the 1-of-K coding scheme.

        Parameters
        ----------
        y : array or sparse matrix of shape [n_samples,] or \
            [n_samples, n_classes]
            Target values. The 2-d matrix should only contain 0 and 1,
            represents multilabel classification. Sparse matrix can be
            CSR, CSC, COO, DOK, or LIL.

        Returns
        -------
        Y : array or CSR matrix of shape [n_samples, n_classes]
            Shape will be [n_samples, 1] for binary problems.
        """
        return self.fit(y).transform(y)

    def transform(self, y):
        """Transform multi-class labels to binary labels

        The output of transform is sometimes referred to by some authors as
        the 1-of-K coding scheme.

        Parameters
        ----------
        y : array or sparse matrix of shape [n_samples,] or \
            [n_samples, n_classes]
            Target values. The 2-d matrix should only contain 0 and 1,
            represents multilabel classification. Sparse matrix can be
            CSR, CSC, COO, DOK, or LIL.

        Returns
        -------
        Y : numpy array or CSR matrix of shape [n_samples, n_classes]
            Shape will be [n_samples, 1] for binary problems.
        """
        check_is_fitted(self, 'classes_')

        y_is_multilabel = type_of_target(y).startswith('multilabel')
        if y_is_multilabel and not self.y_type_.startswith('multilabel'):
            raise ValueError("The object was not fitted with multilabel"
                             " input.")

        return label_binarize(y, self.classes_,
                              pos_label=self.pos_label,
                              neg_label=self.neg_label,
                              sparse_output=self.sparse_output)

    def inverse_transform(self, Y, threshold=None):
        """Transform binary labels back to multi-class labels

        Parameters
        ----------
        Y : numpy array or sparse matrix with shape [n_samples, n_classes]
            Target values. All sparse matrices are converted to CSR before
            inverse transformation.

        threshold : float or None
            Threshold used in the binary and multi-label cases.

            Use 0 when:
                - Y contains the output of decision_function (classifier)
            Use 0.5 when:
                - Y contains the output of predict_proba

            If None, the threshold is assumed to be half way between
            neg_label and pos_label.

        Returns
        -------
        y : numpy array or CSR matrix of shape [n_samples] Target values.

        Notes
        -----
        In the case when the binary labels are fractional
        (probabilistic), inverse_transform chooses the class with the
        greatest value. Typically, this allows to use the output of a
        linear model's decision_function method directly as the input
        of inverse_transform.
        """
        check_is_fitted(self, 'classes_')

        if threshold is None:
            threshold = (self.pos_label + self.neg_label) / 2.

        if self.y_type_ == "multiclass":
            y_inv = _inverse_binarize_multiclass(Y, self.classes_)
        else:
            y_inv = _inverse_binarize_thresholding(Y, self.y_type_,
                                                   self.classes_, threshold)

        if self.sparse_input_:
            y_inv = sp.csr_matrix(y_inv)
        elif sp.issparse(y_inv):
            y_inv = y_inv.toarray()

        return y_inv


def label_binarize(y, classes, neg_label=0, pos_label=1, sparse_output=False):
    """Binarize labels in a one-vs-all fashion

    Several regression and binary classification algorithms are
    available in the scikit. A simple way to extend these algorithms
    to the multi-class classification case is to use the so-called
    one-vs-all scheme.

    This function makes it possible to compute this transformation for a
    fixed set of class labels known ahead of time.

    Parameters
    ----------
    y : array-like
        Sequence of integer labels or multilabel data to encode.

    classes : array-like of shape [n_classes]
        Uniquely holds the label for each class.

    neg_label : int (default: 0)
        Value with which negative labels must be encoded.

    pos_label : int (default: 1)
        Value with which positive labels must be encoded.

    sparse_output : boolean (default: False),
        Set to true if output binary array is desired in CSR sparse format

    Returns
    -------
    Y : numpy array or CSR matrix of shape [n_samples, n_classes]
        Shape will be [n_samples, 1] for binary problems.

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

    See also
    --------
    LabelBinarizer : class used to wrap the functionality of label_binarize and
        allow for fitting to classes independently of the transform operation
    """
    if not isinstance(y, list):
        # XXX Workaround that will be removed when list of list format is
        # dropped
        y = check_array(y, accept_sparse='csr', ensure_2d=False, dtype=None)
    else:
        if _num_samples(y) == 0:
            raise ValueError('y has 0 samples: %r' % y)
    if neg_label >= pos_label:
        raise ValueError("neg_label={0} must be strictly less than "
                         "pos_label={1}.".format(neg_label, pos_label))

    if (sparse_output and (pos_label == 0 or neg_label != 0)):
        raise ValueError("Sparse binarization is only supported with non "
                         "zero pos_label and zero neg_label, got "
                         "pos_label={0} and neg_label={1}"
                         "".format(pos_label, neg_label))

    # To account for pos_label == 0 in the dense case
    pos_switch = pos_label == 0
    if pos_switch:
        pos_label = -neg_label

    y_type = type_of_target(y)
    if 'multioutput' in y_type:
        raise ValueError("Multioutput target data is not supported with label "
                         "binarization")
    if y_type == 'unknown':
        raise ValueError("The type of target data is not known")

    n_samples = y.shape[0] if sp.issparse(y) else len(y)
    n_classes = len(classes)
    classes = np.asarray(classes)

    if y_type == "binary":
        if n_classes == 1:
            if sparse_output:
                return sp.csr_matrix((n_samples, 1), dtype=int)
            else:
                Y = np.zeros((len(y), 1), dtype=np.int)
                Y += neg_label
                return Y
        elif len(classes) >= 3:
            y_type = "multiclass"

    sorted_class = np.sort(classes)
    if (y_type == "multilabel-indicator" and classes.size != y.shape[1]):
        raise ValueError("classes {0} missmatch with the labels {1}"
                         "found in the data".format(classes, unique_labels(y)))

    if y_type in ("binary", "multiclass"):
        y = column_or_1d(y)

        # pick out the known labels from y
        y_in_classes = in1d(y, classes)
        y_seen = y[y_in_classes]
        indices = np.searchsorted(sorted_class, y_seen)
        indptr = np.hstack((0, np.cumsum(y_in_classes)))

        data = np.empty_like(indices)
        data.fill(pos_label)
        Y = sp.csr_matrix((data, indices, indptr),
                          shape=(n_samples, n_classes))
    elif y_type == "multilabel-indicator":
        Y = sp.csr_matrix(y)
        if pos_label != 1:
            data = np.empty_like(Y.data)
            data.fill(pos_label)
            Y.data = data
    else:
        raise ValueError("%s target data is not supported with label "
                         "binarization" % y_type)

    if not sparse_output:
        Y = Y.toarray()
        Y = astype(Y, int, copy=False)

        if neg_label != 0:
            Y[Y == 0] = neg_label

        if pos_switch:
            Y[Y == pos_label] = 0
    else:
        Y.data = astype(Y.data, int, copy=False)

    # preserve label ordering
    if np.any(classes != sorted_class):
        indices = np.searchsorted(sorted_class, classes)
        Y = Y[:, indices]

    if y_type == "binary":
        if sparse_output:
            Y = Y.getcol(-1)
        else:
            Y = Y[:, -1].reshape((-1, 1))

    return Y


def _inverse_binarize_multiclass(y, classes):
    """Inverse label binarization transformation for multiclass.

    Multiclass uses the maximal score instead of a threshold.
    """
    classes = np.asarray(classes)

    if sp.issparse(y):
        # Find the argmax for each row in y where y is a CSR matrix

        y = y.tocsr()
        n_samples, n_outputs = y.shape
        outputs = np.arange(n_outputs)
        row_max = sparse_min_max(y, 1)[1]
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
        samples = np.arange(n_samples)[(row_nnz > 0) &
                                       (row_max.ravel() == 0)]
        for i in samples:
            ind = y.indices[y.indptr[i]:y.indptr[i + 1]]
            y_i_argmax[i] = classes[np.setdiff1d(outputs, ind)][0]

        return classes[y_i_argmax]
    else:
        return classes.take(y.argmax(axis=1), mode="clip")


def _inverse_binarize_thresholding(y, output_type, classes, threshold):
    """Inverse label binarization transformation using thresholding."""

    if output_type == "binary" and y.ndim == 2 and y.shape[1] > 2:
        raise ValueError("output_type='binary', but y.shape = {0}".
                         format(y.shape))

    if output_type != "binary" and y.shape[1] != len(classes):
        raise ValueError("The number of class is not equal to the number of "
                         "dimension of y.")

    classes = np.asarray(classes)

    # Perform thresholding
    if sp.issparse(y):
        if threshold > 0:
            if y.format not in ('csr', 'csc'):
                y = y.tocsr()
            y.data = np.array(y.data > threshold, dtype=np.int)
            y.eliminate_zeros()
        else:
            y = np.array(y.toarray() > threshold, dtype=np.int)
    else:
        y = np.array(y > threshold, dtype=np.int)

    # Inverse transform data
    if output_type == "binary":
        if sp.issparse(y):
            y = y.toarray()
        if y.ndim == 2 and y.shape[1] == 2:
            return classes[y[:, 1]]
        else:
            if len(classes) == 1:
                return np.repeat(classes[0], len(y))
            else:
                return classes[y.ravel()]

    elif output_type == "multilabel-indicator":
        return y

    else:
        raise ValueError("{0} format is not supported".format(output_type))


class MultiLabelBinarizer(BaseEstimator, TransformerMixin):
    """Transform between iterable of iterables and a multilabel format

    Although a list of sets or tuples is a very intuitive format for multilabel
    data, it is unwieldy to process. This transformer converts between this
    intuitive format and the supported multilabel format: a (samples x classes)
    binary matrix indicating the presence of a class label.

    Parameters
    ----------
    classes : array-like of shape [n_classes] (optional)
        Indicates an ordering for the class labels

    sparse_output : boolean (default: False),
        Set to true if output binary array is desired in CSR sparse format

    Attributes
    ----------
    classes_ : array of labels
        A copy of the `classes` parameter where provided,
        or otherwise, the sorted set of classes found when fitting.

    Examples
    --------
    >>> from sklearn.preprocessing import MultiLabelBinarizer
    >>> mlb = MultiLabelBinarizer()
    >>> mlb.fit_transform([(1, 2), (3,)])
    array([[1, 1, 0],
           [0, 0, 1]])
    >>> mlb.classes_
    array([1, 2, 3])

    >>> mlb.fit_transform([set(['sci-fi', 'thriller']), set(['comedy'])])
    array([[0, 1, 1],
           [1, 0, 0]])
    >>> list(mlb.classes_)
    ['comedy', 'sci-fi', 'thriller']

    See also
    --------
    sklearn.preprocessing.OneHotEncoder : encode categorical integer features
        using a one-hot aka one-of-K scheme.
    """
    def __init__(self, classes=None, sparse_output=False):
        self.classes = classes
        self.sparse_output = sparse_output

    def fit(self, y):
        """Fit the label sets binarizer, storing `classes_`

        Parameters
        ----------
        y : iterable of iterables
            A set of labels (any orderable and hashable object) for each
            sample. If the `classes` parameter is set, `y` will not be
            iterated.

        Returns
        -------
        self : returns this MultiLabelBinarizer instance
        """
        if self.classes is None:
            classes = sorted(set(itertools.chain.from_iterable(y)))
        else:
            classes = self.classes
        dtype = np.int if all(isinstance(c, int) for c in classes) else object
        self.classes_ = np.empty(len(classes), dtype=dtype)
        self.classes_[:] = classes
        return self

    def fit_transform(self, y):
        """Fit the label sets binarizer and transform the given label sets

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
        if self.classes is not None:
            return self.fit(y).transform(y)

        # Automatically increment on new class
        class_mapping = defaultdict(int)
        class_mapping.default_factory = class_mapping.__len__
        yt = self._transform(y, class_mapping)

        # sort classes and reorder columns
        tmp = sorted(class_mapping, key=class_mapping.get)

        # (make safe for tuples)
        dtype = np.int if all(isinstance(c, int) for c in tmp) else object
        class_mapping = np.empty(len(tmp), dtype=dtype)
        class_mapping[:] = tmp
        self.classes_, inverse = np.unique(class_mapping, return_inverse=True)
        # ensure yt.indices keeps its current dtype
        yt.indices = np.array(inverse[yt.indices], dtype=yt.indices.dtype,
                              copy=False)

        if not self.sparse_output:
            yt = yt.toarray()

        return yt

    def transform(self, y):
        """Transform the given label sets

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
        check_is_fitted(self, 'classes_')

        class_to_index = dict(zip(self.classes_, range(len(self.classes_))))
        yt = self._transform(y, class_to_index)

        if not self.sparse_output:
            yt = yt.toarray()

        return yt

    def _transform(self, y, class_mapping):
        """Transforms the label sets with a given mapping

        Parameters
        ----------
        y : iterable of iterables
        class_mapping : Mapping
            Maps from label to column index in label indicator matrix

        Returns
        -------
        y_indicator : sparse CSR matrix, shape (n_samples, n_classes)
            Label indicator matrix
        """
        indices = array.array('i')
        indptr = array.array('i', [0])
        for labels in y:
            indices.extend(set(class_mapping[label] for label in labels))
            indptr.append(len(indices))
        data = np.ones(len(indices), dtype=int)

        return sp.csr_matrix((data, indices, indptr),
                             shape=(len(indptr) - 1, len(class_mapping)))

    def inverse_transform(self, yt):
        """Transform the given indicator matrix into label sets

        Parameters
        ----------
        yt : array or sparse matrix of shape (n_samples, n_classes)
            A matrix containing only 1s ands 0s.

        Returns
        -------
        y : list of tuples
            The set of labels for each sample such that `y[i]` consists of
            `classes_[j]` for each `yt[i, j] == 1`.
        """
        check_is_fitted(self, 'classes_')

        if yt.shape[1] != len(self.classes_):
            raise ValueError('Expected indicator for {0} classes, but got {1}'
                             .format(len(self.classes_), yt.shape[1]))

        if sp.issparse(yt):
            yt = yt.tocsr()
            if len(yt.data) != 0 and len(np.setdiff1d(yt.data, [0, 1])) > 0:
                raise ValueError('Expected only 0s and 1s in label indicator.')
            return [tuple(self.classes_.take(yt.indices[start:end]))
                    for start, end in zip(yt.indptr[:-1], yt.indptr[1:])]
        else:
            unexpected = np.setdiff1d(yt, [0, 1])
            if len(unexpected) > 0:
                raise ValueError('Expected only 0s and 1s in label indicator. '
                                 'Also got {0}'.format(unexpected))
            return [tuple(self.classes_.compress(indicators)) for indicators
                    in yt]
