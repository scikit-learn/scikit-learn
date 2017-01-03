from ..base import BaseEstimator, TransformerMixin
from ..utils import check_array


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

    kw_args : dict, optional
        Dictionary of additional keyword arguments to pass to func.

    inv_kw_args : dict, optional
        Dictionary of additional keyword arguments to pass to inverse_func.

    """
    def __init__(self, func=None, inverse_func=None, validate=True,
                 accept_sparse=False, pass_y=False,
                 kw_args=None, inv_kw_args=None):
        self.func = func
        self.inverse_func = inverse_func
        self.validate = validate
        self.accept_sparse = accept_sparse
        self.pass_y = pass_y
        self.kw_args = kw_args
        self.inv_kw_args = inv_kw_args

    def fit(self, X, y=None):
        if self.validate:
            check_array(X, self.accept_sparse)
        return self

    def transform(self, X, y=None):
        return self._transform(X, y, self.func, self.kw_args)

    def inverse_transform(self, X, y=None):
        return self._transform(X, y, self.inverse_func, self.inv_kw_args)

    def _transform(self, X, y=None, func=None, kw_args=None):
        if self.validate:
            X = check_array(X, self.accept_sparse)

        if func is None:
            func = _identity

        return func(X, *((y,) if self.pass_y else ()),
                    **(kw_args if kw_args else {}))
