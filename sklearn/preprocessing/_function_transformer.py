import warnings

from ..base import BaseEstimator, TransformerMixin
from ..utils.validation import (
    _allclose_dense_sparse,
    _check_feature_names_in,
    check_array,
    column_or_1d,
)


def _identity(X):
    """The identity function."""
    return X


class FunctionTransformer(TransformerMixin, BaseEstimator):
    """Constructs a transformer from an arbitrary callable.

    A FunctionTransformer forwards its X (and optionally y) arguments to a
    user-defined function or function object and returns the result of this
    function. This is useful for stateless transformations such as taking the
    log of frequencies, doing custom scaling, etc.

    Note: If a lambda is used as the function, then the resulting
    transformer will not be pickleable.

    .. versionadded:: 0.17

    Read more in the :ref:`User Guide <function_transformer>`.

    Parameters
    ----------
    func : callable, default=None
        The callable to use for the transformation. This will be passed
        the same arguments as transform, with args and kwargs forwarded.
        If func is None, then func will be the identity function.

    inverse_func : callable, default=None
        The callable to use for the inverse transformation. This will be
        passed the same arguments as inverse transform, with args and
        kwargs forwarded. If inverse_func is None, then inverse_func
        will be the identity function.

    validate : bool, default=False
        Indicate that the input X array should be checked before calling
        ``func``. The possibilities are:

        - If False, there is no input validation.
        - If True, then X will be converted to a 2-dimensional NumPy array or
          sparse matrix. If the conversion is not possible an exception is
          raised.

        .. versionchanged:: 0.22
           The default of ``validate`` changed from True to False.

    accept_sparse : bool, default=False
        Indicate that func accepts a sparse matrix as input. If validate is
        False, this has no effect. Otherwise, if accept_sparse is false,
        sparse matrix inputs will cause an exception to be raised.

    check_inverse : bool, default=True
       Whether to check that or ``func`` followed by ``inverse_func`` leads to
       the original inputs. It can be used for a sanity check, raising a
       warning when the condition is not fulfilled.

       .. versionadded:: 0.20

    feature_names_out : array-like of str, or callable, or None, default=None
        Determines the list of feature names that will be returned by the
        `get_feature_names_out` method. If you pass a callable, then it must
        take two positional arguments: this `FunctionTransformer` (`self`) and
        an array-like of input feature names (`input_features`). It must return
        an array-like of output feature names.

        See ``get_feature_names_out`` for more details.

        .. versionadded:: 1.1

    kw_args : dict, default=None
        Dictionary of additional keyword arguments to pass to func.

        .. versionadded:: 0.18

    inv_kw_args : dict, default=None
        Dictionary of additional keyword arguments to pass to inverse_func.

        .. versionadded:: 0.18

    Attributes
    ----------
    n_features_in_ : int
        Number of features seen during :term:`fit`. Defined only when
        `validate=True`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `validate=True`
        and `X` has feature names that are all strings.

        .. versionadded:: 1.0

    See Also
    --------
    MaxAbsScaler : Scale each feature by its maximum absolute value.
    StandardScaler : Standardize features by removing the mean and
        scaling to unit variance.
    LabelBinarizer : Binarize labels in a one-vs-all fashion.
    MultiLabelBinarizer : Transform between iterable of iterables
        and a multilabel format.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.preprocessing import FunctionTransformer
    >>> transformer = FunctionTransformer(np.log1p)
    >>> X = np.array([[0, 1], [2, 3]])
    >>> transformer.transform(X)
    array([[0.       , 0.6931...],
           [1.0986..., 1.3862...]])
    """

    def __init__(
        self,
        func=None,
        inverse_func=None,
        *,
        validate=False,
        accept_sparse=False,
        check_inverse=True,
        feature_names_out=None,
        kw_args=None,
        inv_kw_args=None,
    ):
        self.func = func
        self.inverse_func = inverse_func
        self.validate = validate
        self.accept_sparse = accept_sparse
        self.check_inverse = check_inverse
        self.feature_names_out = feature_names_out
        self.kw_args = kw_args
        self.inv_kw_args = inv_kw_args

    def _check_input(self, X, *, reset):
        if self.validate:
            return self._validate_data(X, accept_sparse=self.accept_sparse, reset=reset)
        return X

    def _check_inverse_transform(self, X):
        """Check that func and inverse_func are the inverse."""
        idx_selected = slice(None, None, max(1, X.shape[0] // 100))
        X_round_trip = self.inverse_transform(self.transform(X[idx_selected]))
        if not _allclose_dense_sparse(X[idx_selected], X_round_trip):
            warnings.warn(
                "The provided functions are not strictly"
                " inverse of each other. If you are sure you"
                " want to proceed regardless, set"
                " 'check_inverse=False'.",
                UserWarning,
            )

    def fit(self, X, y=None):
        """Fit transformer by checking X.

        If ``validate`` is ``True``, ``X`` will be checked.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Input array.

        y : Ignored
            Not used, present here for API consistency by convention.

        Returns
        -------
        self : object
            FunctionTransformer class instance.
        """
        X = self._check_input(X, reset=True)
        if self.check_inverse and not (self.func is None or self.inverse_func is None):
            self._check_inverse_transform(X)
        return self

    def transform(self, X):
        """Transform X using the forward function.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Input array.

        Returns
        -------
        X_out : array-like, shape (n_samples, n_features)
            Transformed input.
        """
        X = self._check_input(X, reset=False)
        return self._transform(X, func=self.func, kw_args=self.kw_args)

    def inverse_transform(self, X):
        """Transform X using the inverse function.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Input array.

        Returns
        -------
        X_out : array-like, shape (n_samples, n_features)
            Transformed input.
        """
        if self.validate:
            X = check_array(X, accept_sparse=self.accept_sparse)
        return self._transform(X, func=self.inverse_func, kw_args=self.inv_kw_args)

    def get_feature_names_out(self, input_features=None):
        """Get output feature names for transformation.

        Parameters
        ----------
        input_features : array-like of str or None, default=None
            Input feature names.

            - If `input_features` is None, then `feature_names_in_` is
              used as the input feature names. If `feature_names_in_` is not
              defined, then names are generated:
              `[x0, x1, ..., x(n_features_in_)]`.
            - If `input_features` is array-like, then `input_features` must
              match `feature_names_in_` if `feature_names_in_` is defined.

        Returns
        -------
        feature_names_out : ndarray of str objects
            Transformed feature names.

            - If `feature_names_out` is None, the input feature names are
              returned (see `input_features` above). This requires
              `n_features_in_` to be defined, which in turn requires
              `validate=True`.
            - If `feature_names_out` is an array-like of strings, then it
              is returned, ignoring the input feature names.
            - If `feature_names_out` is a callable, then it is called with two
              arguments, `self` and `input_features`, and its return value is
              returned by this method.
        """
        if hasattr(self, "n_features_in_") or input_features is not None:
            input_features = _check_feature_names_in(self, input_features)
        elif input_features is not None:
            input_features = column_or_1d(input_features)
        if self.feature_names_out is None:
            if input_features is None:
                raise ValueError(
                    "If 'feature_names_out' is None, then 'input_features' "
                    "must be passed, or 'n_features_in_' must be defined. If "
                    "you set 'validate' to 'True', then 'n_features_in_' will "
                    "be set automatically when 'fit' is called."
                )
            names_out = input_features
        elif callable(self.feature_names_out):
            names_out = self.feature_names_out(self, input_features)
        elif isinstance(self.feature_names_out, str):
            raise ValueError(
                "'feature_names_out' must not be a string. If there is a "
                "single output feature name, then set 'feature_names_out' to "
                "an array-like containing just that name."
            )
        else:
            names_out = self.feature_names_out
        return column_or_1d(names_out)

    def _transform(self, X, func=None, kw_args=None):
        if func is None:
            func = _identity

        return func(X, **(kw_args if kw_args else {}))

    def __sklearn_is_fitted__(self):
        """Return True since FunctionTransfomer is stateless."""
        return True

    def _more_tags(self):
        return {"no_validation": not self.validate, "stateless": True}
