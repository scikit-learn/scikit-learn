import warnings

from ..base import BaseEstimator, TransformerMixin
from ..utils import check_array
from ..utils.validation import _allclose_dense_sparse


def _identity(X):
    """The identity function.
    """
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
    func : callable, optional default=None
        The callable to use for the transformation. This will be passed
        the same arguments as transform, with args and kwargs forwarded.
        If func is None, then func will be the identity function.

    inverse_func : callable, optional default=None
        The callable to use for the inverse transformation. This will be
        passed the same arguments as inverse transform, with args and
        kwargs forwarded. If inverse_func is None, then inverse_func
        will be the identity function.

    validate : bool, optional default=False
        Indicate that the input X array should be checked before calling
        ``func``. The possibilities are:

        - If False, there is no input validation.
        - If True, then X will be converted to a 2-dimensional NumPy array or
          sparse matrix. If the conversion is not possible an exception is
          raised.

        .. deprecated:: 0.22
           The default of ``validate`` changed from True to False.

    accept_sparse : boolean, optional
        Indicate that func accepts a sparse matrix as input. If validate is
        False, this has no effect. Otherwise, if accept_sparse is false,
        sparse matrix inputs will cause an exception to be raised.

    check_inverse : bool, default=True
       Whether to check that or ``func`` followed by ``inverse_func`` leads to
       the original inputs. It can be used for a sanity check, raising a
       warning when the condition is not fulfilled.

       .. versionadded:: 0.20

    kw_args : dict, optional
        Dictionary of additional keyword arguments to pass to func.

    inv_kw_args : dict, optional
        Dictionary of additional keyword arguments to pass to inverse_func.

    """
    def __init__(self, func=None, inverse_func=None, validate=False,
                 accept_sparse=False, check_inverse=True, kw_args=None,
                 inv_kw_args=None):
        self.func = func
        self.inverse_func = inverse_func
        self.validate = validate
        self.accept_sparse = accept_sparse
        self.check_inverse = check_inverse
        self.kw_args = kw_args
        self.inv_kw_args = inv_kw_args

    def _check_input(self, X):
        if self.validate:
            return check_array(X, accept_sparse=self.accept_sparse)
        return X

    def _check_inverse_transform(self, X):
        """Check that func and inverse_func are the inverse."""
        idx_selected = slice(None, None, max(1, X.shape[0] // 100))
        X_round_trip = self.inverse_transform(self.transform(X[idx_selected]))
        if not _allclose_dense_sparse(X[idx_selected], X_round_trip):
            warnings.warn("The provided functions are not strictly"
                          " inverse of each other. If you are sure you"
                          " want to proceed regardless, set"
                          " 'check_inverse=False'.", UserWarning)

    def fit(self, X, y=None):
        """Fit transformer by checking X.

        If ``validate`` is ``True``, ``X`` will be checked.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Input array.

        Returns
        -------
        self
        """
        X = self._check_input(X)
        if (self.check_inverse and not (self.func is None or
                                        self.inverse_func is None)):
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
        return self._transform(X, func=self.inverse_func,
                               kw_args=self.inv_kw_args)

    def _transform(self, X, func=None, kw_args=None):
        X = self._check_input(X)

        if func is None:
            func = _identity

        return func(X, **(kw_args if kw_args else {}))

    def _more_tags(self):
        return {'no_validation': not self.validate,
                'stateless': True}
