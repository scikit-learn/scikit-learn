import warnings

from ..base import BaseEstimator, TransformerMixin
from ..utils import check_array
from ..utils.testing import assert_allclose_dense_sparse
from ..utils.validation import _assert_all_finite
from ..externals.six import string_types


def _identity(X):
    """The identity function.
    """
    return X


class FunctionTransformer(BaseEstimator, TransformerMixin):
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

    validate : bool or 'array-or-frame', optional default=True
        Indicate that the input X array should be checked before calling
        func. The possibilities are:

        - If 'array-or-frame', X will be passed through if it is a pandas
          DataFrame or converted to a 2-dimensional array or sparse matrix. In
          this latest case, an exception will be raised if the conversion
          failed.
        - If True, then X will be converted to a 2-dimensional NumPy array or
          sparse matrix. If the conversion is not possible an exception is
          raised.
        - If False, there is no input validation.

        When X is validated, the parameters ``accept_sparse`` and
        ``force_all_finite`` control the validation for the sparsity and
        the finiteness of X, respectively.

        .. deprecated:: 0.20
           ``validate=True`` as default will be replaced by
           ``validate='array-or-frame'`` in 0.22.

        .. versionchanged:: 0.20
           ``validate`` takes the option ``'array-or-frame'``.

    accept_sparse : boolean, optional
        Indicate that func accepts a sparse matrix as input. If validate is
        False, this has no effect. Otherwise, if accept_sparse is false,
        sparse matrix inputs will cause an exception to be raised.

    force_all_finite : boolean or 'allow-nan', optional default=True
        Whether to raise an error on np.inf and np.nan in X. If validate is
        False, this has not effect. The possibilities are:

        - If True, force all values of X to be finite.
        - If False, accept both np.inf and np.nan in X.
        - If 'allow-nan', accept only np.nan values in X. Values cannot be
          infinite.

        .. versionadded:: 0.20

    pass_y : bool, optional default=False
        Indicate that transform should forward the y argument to the
        inner callable.

        .. deprecated::0.19

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
    def __init__(self, func=None, inverse_func=None, validate=None,
                 accept_sparse=False, force_all_finite=True,
                 pass_y='deprecated', check_inverse=True, kw_args=None,
                 inv_kw_args=None):
        self.func = func
        self.inverse_func = inverse_func
        self.validate = validate
        self.accept_sparse = accept_sparse
        self.force_all_finite = force_all_finite
        self.pass_y = pass_y
        self.check_inverse = check_inverse
        self.kw_args = kw_args
        self.inv_kw_args = inv_kw_args

    def _check_input(self, X):
        # FIXME: Future warning to be removed in 0.22
        if self.validate is None:
            self._validate = True
            if hasattr(X, 'loc'):
                warnings.warn("The default validate=True will be replaced by "
                              "validate='array-or-frame' in 0.22. A pandas "
                              "DataFrame will not be converted to a 2D "
                              "NumPy array.", FutureWarning)
        else:
            self._validate = self.validate

        if self._validate not in (True, False, 'array-or-frame'):
            raise ValueError("'validate' should be a boolean or "
                             "'array-or-frame'. Got {!r} instead."
                             .format(self._validate))
        if ((not isinstance(self.force_all_finite, bool)) and
                self.force_all_finite != 'allow-nan'):
            raise ValueError("'force_all_finite' should be a boolean "
                             "or 'allow-nan'. Got {!r} instead."
                             .format(self.force_all_finite))

        if self._validate:
            if hasattr(X, 'loc') and self._validate == 'array-or-frame':
                if self.force_all_finite:
                    _assert_all_finite(
                        X.values,
                        allow_nan=not (self.force_all_finite is True))
                return X
            else:
                return check_array(X, accept_sparse=self.accept_sparse,
                                   force_all_finite=self.force_all_finite)
        return X

    def _check_inverse_transform(self, X):
        """Check that func and inverse_func are the inverse."""
        idx_selected = slice(None, None, max(1, X.shape[0] // 100))
        try:
            assert_allclose_dense_sparse(
                X[idx_selected],
                self.inverse_transform(self.transform(X[idx_selected])))
        except AssertionError:
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

    def transform(self, X, y='deprecated'):
        """Transform X using the forward function.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Input array.

        y : (ignored)
            .. deprecated::0.19

        Returns
        -------
        X_out : array-like, shape (n_samples, n_features)
            Transformed input.
        """
        if not isinstance(y, string_types) or y != 'deprecated':
            warnings.warn("The parameter y on transform() is "
                          "deprecated since 0.19 and will be removed in 0.21",
                          DeprecationWarning)

        return self._transform(X, y=y, func=self.func, kw_args=self.kw_args)

    def inverse_transform(self, X, y='deprecated'):
        """Transform X using the inverse function.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Input array.

        y : (ignored)
            .. deprecated::0.19

        Returns
        -------
        X_out : array-like, shape (n_samples, n_features)
            Transformed input.
        """
        if not isinstance(y, string_types) or y != 'deprecated':
            warnings.warn("The parameter y on inverse_transform() is "
                          "deprecated since 0.19 and will be removed in 0.21",
                          DeprecationWarning)
        return self._transform(X, y=y, func=self.inverse_func,
                               kw_args=self.inv_kw_args)

    def _transform(self, X, y=None, func=None, kw_args=None):
        X = self._check_input(X)

        if func is None:
            func = _identity

        if (not isinstance(self.pass_y, string_types) or
                self.pass_y != 'deprecated'):
            # We do this to know if pass_y was set to False / True
            pass_y = self.pass_y
            warnings.warn("The parameter pass_y is deprecated since 0.19 and "
                          "will be removed in 0.21", DeprecationWarning)
        else:
            pass_y = False

        return func(X, *((y,) if pass_y else ()),
                    **(kw_args if kw_args else {}))
