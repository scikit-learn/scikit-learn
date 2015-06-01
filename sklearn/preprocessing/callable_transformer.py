from ..base import BaseEstimator, TransformerMixin
from ..utils import check_array


class CallableTransformer(BaseEstimator, TransformerMixin):
    """Allows the construction of a transformer from an arbitrary callable.

    Parameters
    ----------
    func : callable, optional default=None
        The callable to use for the transformation. This will be passed
        the same arguments as transform, with args and kwargs forwarded.
        If func is None, then func will be the identity function.
    validate : bool, optional default=True
        Indicate that the input X array should be checked before calling
        func. If validate is false, there will be no input validation.
    accept_sparse : boolean, optional
        Indicate that func accepts a sparse matrix as input.
    args : tuple, optional
        A tuple of positional arguments to be passed to func. These will
        be passed after X and y.
    kwargs : dict, optional
        A dictionary of keyword arguments to be passed to func.

    """
    def __init__(self, func=None, validate=True, accept_sparse=False,
                 args=None, kwargs=None):
        self.func = func
        self.validate = validate
        self.accept_sparse = accept_sparse
        self.args = args
        self.kwargs = kwargs

    def fit(self, X, y=None):
        if self.validate:
            check_array(X, self.accept_sparse)
        return self

    def transform(self, X, y=None):
        if self.validate:
            X = check_array(X, self.accept_sparse)
        return (self.func or (lambda X, y, *args, **kwargs: X))(
            X, y, *(self.args or ()), **(self.kwargs or {})
        )
