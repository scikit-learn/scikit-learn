import warnings

import numpy as np

from ..base import BaseEstimator, TransformerMixin
from ..utils import check_array, check_random_state
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

    A FunctionTransformer will not do any checks on its function's output.

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

    validate : bool, optional default=True
        Indicate that the input X array should be checked before calling
        func. If validate is false, there will be no input validation.
        If it is true, then X will be converted to a 2-dimensional NumPy
        array or sparse matrix. If this conversion is not possible or X
        contains NaN or infinity, an exception is raised.

    accept_sparse : boolean, optional
        Indicate that func accepts a sparse matrix as input. If validate is
        False, this has no effect. Otherwise, if accept_sparse is false,
        sparse matrix inputs will cause an exception to be raised.

    pass_y : bool, optional default=False
        Indicate that transform should forward the y argument to the
        inner callable.

        .. deprecated::0.19

    check_inverse : bool, (default=False)
       Whether to check that ``transform`` followed by ``inverse_transform``
       or ``func`` followed by ``inverse_func`` leads to the original targets.

       .. versionadded:: 0.20

    kw_args : dict, optional
        Dictionary of additional keyword arguments to pass to func.

    inv_kw_args : dict, optional
        Dictionary of additional keyword arguments to pass to inverse_func.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by np.random. Note that this is used to compute if func and
        inverse_func are the inverse of each other.


    """
    def __init__(self, func=None, inverse_func=None, validate=True,
                 accept_sparse=False, pass_y='deprecated', check_inverse=False,
                 kw_args=None, inv_kw_args=None, random_state=None):
        self.func = func
        self.inverse_func = inverse_func
        self.validate = validate
        self.accept_sparse = accept_sparse
        self.pass_y = pass_y
        self.check_inverse = check_inverse
        self.kw_args = kw_args
        self.inv_kw_args = inv_kw_args
        self.random_state = random_state

    def _validate_inverse(self, X):
        """Check that func and inverse_func are the inverse."""
        # Apply the transform and inverse_transform on few samples.
        random_state = check_random_state(self.random_state)
        n_subsample = min(10, X.shape[0])
        subsample_idx = random_state.choice(range(X.shape[0]),
                                            size=n_subsample,
                                            replace=False)
        if not np.allclose(X[subsample_idx],
                           self.inverse_transform(
                               self.transform(X[subsample_idx])),
                           atol=1e-7):
            raise ValueError("The provided functions are not strictly"
                             " inverse of each other. If you are sure you"
                             " want to proceed regardless, set"
                             " 'check_inverse=False'")

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
        if self.validate:
            check_array(X, self.accept_sparse)
        if self.check_inverse:
            self._validate_inverse(X)
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
        if self.validate:
            X = check_array(X, self.accept_sparse)

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
