from functools import wraps

from scipy.sparse import issparse

from sklearn.utils import check_pandas_support
from sklearn._config import get_config

__all__ = [
    "get_output_config",
    "SetOutputMixin",
    "safe_set_output",
]


def _wrap_in_pandas_container(
    original_data,
    *,
    index=None,
    columns=None,
    constructor_kwargs=None,
):
    """Create a named container.

    Parameters
    ----------
    original_data : ndarray, sparse matrix or pandas DataFrame
        Container to name.

    index : array-like, default=None
        Index for data.

    columns : callable or ndarray, default=None
        The column names or a callable that returns the column names. This is
        useful if the column names require some computation.

    constructor_kwargs : dict, default=None
        Keyword arguments passed to container constructor.

    Returns
    -------
    named_container : DataFrame or ndarray
        Container with column names or unchanged `output`.
    """
    # Already a pandas DataFrame
    if hasattr(original_data, "iloc"):
        return original_data

    constructor_kwargs = constructor_kwargs or {}

    if issparse(original_data):
        raise ValueError("Pandas output does not support sparse data")

    if callable(columns):
        columns = columns()

    pd = check_pandas_support("Setting output container to 'pandas'")
    return pd.DataFrame(
        original_data, index=index, columns=columns, **constructor_kwargs
    )


def get_output_config(estimator, method):
    """Get output configure based on estimator and global configuration.

    Parameters
    ----------
    estimator : estimator instance
        If not `None`, check the estimator for output container.

    method : {"transform"}
        Method to get container output for.

    Returns
    -------
    config : dict
        Dictionary with keys, "dense", that specifies the container for `method`.
    """
    est_sklearn_output_config = getattr(estimator, "_sklearn_output_config", {})
    if method in est_sklearn_output_config:
        container_str = est_sklearn_output_config[method]
    else:
        container_str = get_config()[f"output_{method}"]

    return {"dense": container_str}


def _wrap_output_with_container(
    estimator, original_data, method, index, constructor_kwargs=None
):
    """Wrap output with container based on an estimator's or global config.

    Parameters
    ----------
    estimator : estimator instance
        Estimator to get the output configuration from.

    original_data : ndarray
        Data to wrap with container.

    method : {"transform"}
        Method to get container output for.

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
    output_config = get_output_config(estimator, method)
    dense_config = output_config["dense"]
    if dense_config not in {"default", "pandas"}:
        raise ValueError(
            f"output config must be 'default' or 'pandas' got {dense_config}"
        )

    if dense_config == "default":
        return original_data

    # dense_config == "pandas"
    return _wrap_in_pandas_container(
        original_data=original_data,
        index=index,
        columns=getattr(estimator, "get_feature_names_out", None),
        constructor_kwargs=constructor_kwargs,
    )


def _wrap_method_output(f, method):
    """Wrapper used by SetOutputMixin to automatically wrap methods."""

    @wraps(f)
    def wrapped(self, X, *args, **kwargs):
        original_data = f(self, X, *args, **kwargs)
        return _wrap_output_with_container(
            self, original_data, method, getattr(X, "index", None)
        )

    return wrapped


class SetOutputMixin:
    """Mixin that dynamically wraps methods to return container based on config.

    Currently `SetOutputMixin` wraps `transform` and `fit_transform` and configures
    it based on `set_output` of the global configuration.
    """

    def __init_subclass__(cls, **kwargs):
        # Dynamically wraps `transform` and `fit_transform` and configure it's
        # output based on `set_output`.
        if hasattr(cls, "transform"):
            cls.transform = _wrap_method_output(cls.transform, "transform")
        if hasattr(cls, "fit_transform"):
            cls.fit_transform = _wrap_method_output(cls.fit_transform, "transform")

    def set_output(self, *, transform=None):
        """Set output container.

        .. note:: Experimental API
            The `set_output` API is experimental and subject to change without
            deprecation.

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


def safe_set_output(estimator, *, transform=None):
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
