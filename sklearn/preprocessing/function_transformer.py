from ..base import BaseEstimator, TransformerMixin
from ..utils import check_array


def _identity(X):
    """The identity function.
    """
    return X


class FunctionTransformer(BaseEstimator, TransformerMixin):
    """Constructs a transformer from an arbitrary callable.

    Note: If a lambda is used as the function, then the resulting
    transformer will not be pickleable.

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

    pass_y: bool, optional default=False
        Indicate that transform should forward the y argument to the
        inner callable.

    """
    def __init__(self, func=None, validate=True,
                 accept_sparse=False, pass_y=False):
        self.func = func
        self.validate = validate
        self.accept_sparse = accept_sparse
        self.pass_y = pass_y

    def fit(self, X, y=None):
        if self.validate:
            check_array(X, self.accept_sparse)
        return self

    def transform(self, X, y=None):
        if self.validate:
            X = check_array(X, self.accept_sparse)

        return (self.func or _identity)(X, *((y,) if self.pass_y else ()))
