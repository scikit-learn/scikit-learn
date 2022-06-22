from functools import wraps

from scipy.sparse import issparse

from sklearn.utils import check_pandas_support
from sklearn._config import get_config
from sklearn.utils import Bunch

__all__ = [
    "get_output_config",
    "make_named_container",
    "OutputTypeMixin",
    "safe_set_output",
]


def make_named_container(
    output,
    *,
    index=None,
    columns=None,
    dense_container="default",
    constructor_kwargs=None,
):
    """Create a named container.

    Parameters
    ----------
    output : ndarray, sparse matrix or pandas DataFrame
        Container to name.

    index : array-like, default=None
        Index for data.

    columns : callable or ndarray, default=None
        The column names or a callable that returns the column names. This is
        useful if the column names require some computation.

    dense_container : {"default", "pandas"}, default="default"
        Container used for dense data.

    constructor_kwargs : dict, default=None
        Keyword arguments passed to container constructor.

    Returns
    -------
    output : DataFrame or ndarray
        Container with column names or unchanged `output`.
    """
    if dense_container not in {"default", "pandas"}:
        raise ValueError(
            f"dense_container must be 'default' or 'pandas' got {dense_container}"
        )

    if dense_container == "default":
        return output

    constructor_kwargs = constructor_kwargs or {}

    # dense_container == "pandas"
    if issparse(output):
        raise ValueError("Sparse data does not support pandas output")

    if callable(columns):
        columns = columns()

    if hasattr(output, "columns"):
        if columns is not None:
            output.columns = columns
        if index is not None:
            output.index = index
        return output

    pd = check_pandas_support("make_named_container")
    return pd.DataFrame(output, index=index, columns=columns, **constructor_kwargs)


def get_output_config(estimator, method):
    """Get output configure based on estimator and global configuration.

    Parameters
    ----------
    estimator : estimator instance, default=None
        If not `None`, check the estimator for output container.

    method : str
        Method to get container output for.

    Returns
    -------
    config : Bunch
        Dictionary with keys, "dense", that specifies the
        container for `method`.
    """
    est_sklearn_output_config = getattr(estimator, "_sklearn_output_config", {})
    if method in est_sklearn_output_config:
        container_str = est_sklearn_output_config[method]
    else:
        container_str = get_config()[f"output_{method}"]

    return Bunch(**{"dense": container_str})


def _wrap_output_with_container(
    estimator, output, method, index, constructor_kwargs=None
):
    """Wrap output with container based on an estimator's or global config.

    Parameters
    ----------
    estimator : estimator instance
        Estimator to get the output configuration from.

    output : ndarray
        Output to to wrap with container.

    method : str
        Method that returned `output`.

    index : array-like
        Index to attach to output.

    constructor_kwargs : dict, default=None
        Keyword arguments passed to container constructor.

    Returns
    -------
    wrapped_output : ndarray or DataFrame
        Wrapped output with column names and index or `output` itself if wrapping is
        not configured.
    """
    output_container = get_output_config(estimator, method)
    return make_named_container(
        output=output,
        index=index,
        columns=getattr(estimator, "get_feature_names_out", None),
        dense_container=output_container.dense,
        constructor_kwargs=constructor_kwargs,
    )


def _wrap_method_output(f, method):
    """Wrapper used by OutputTypeMixin to automatically wrap methods."""

    @wraps(f)
    def wrapped(self, X, *args, **kwargs):
        output = f(self, X, *args, **kwargs)
        return _wrap_output_with_container(
            self, output, method, getattr(X, "index", None)
        )

    return wrapped


class OutputTypeMixin:
    """Mixin that dynamically wraps methods to return container based on config.

    Currently `OutputTypeMixin` wraps `transform` and `fit_transform` and configures
    it based on `set_output` of the global configuration.
    """

    def __init_subclass__(cls, **kwargs):
        # Dynamically wraps `transform` and `fit_transform` and configure it's
        # output based on `set_output`.
        if hasattr(cls, "transform"):
            cls.transform = _wrap_method_output(cls.transform, "transform")
        if hasattr(cls, "fit_transform"):
            cls.fit_transform = _wrap_method_output(cls.fit_transform, "transform")

    def set_output(self, transform=None):
        """Set output container.

        Parameters
        ----------
        transform : {"default", "pandas"}, default=None
            Configure output of `transform` and `fit_transform`.

        Returns
        -------
        self : estimator instance
            Estimator instance.
        """
        if transform is None:
            return

        if not hasattr(self, "_sklearn_output_config"):
            self._sklearn_output_config = {}

        self._sklearn_output_config["transform"] = transform
        return self


def safe_set_output(estimator, transform=None):
    """Safely call estimator.set_output and error if it not available.

    This is used by meta-estimators to set the output for child estimators.

    Parameters
    ----------
    estimator : estimator instance
        Estimator instance.

    transform : {"default", "pandas"}, default=None
        Configure output of `transform` and `fit_transform`.

    Returns
    -------
    estimator : estimator instance
        Estimator instance.
    """
    set_output_for_transform = (
        hasattr(estimator, "transform")
        or hasattr(estimator, "fit_transform")
        and transform is not None
    )
    if not set_output_for_transform:
        # If estimator can not transform, then `set_output` does not need to be
        # called.
        return

    if not hasattr(estimator, "set_output"):
        raise ValueError(
            f"Unable to configure output for {estimator} because `set_output` "
            "is not available"
        )
    return estimator.set_output(transform=transform)
