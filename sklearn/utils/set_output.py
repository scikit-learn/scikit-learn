from functools import wraps

from scipy.sparse import issparse

from ..utils import check_pandas_support
from .._config import get_config
from ._available_if import available_if


__all__ = [
    "get_output_config",
    "SetOutputMixin",
    "safe_set_output",
]


def _wrap_in_pandas_container(
    data_to_wrap,
    *,
    index=None,
    columns=None,
):
    """Create a Pandas DataFrame.

    Parameters
    ----------
    data_to_wrap : {ndarray, dataframe}
        Container to name.

    index : array-like, default=None
        Index for data.

    columns : callable or ndarray, default=None
        The column names or a callable that returns the column names. This is
        useful if the column names require some computation.

    Returns
    -------
    dataframe : DataFrame
        Container with column names or unchanged `output`.
    """
    if issparse(data_to_wrap):
        raise ValueError("Pandas output does not support sparse data")

    if callable(columns):
        columns = columns()

    # Already a pandas DataFrame
    if hasattr(data_to_wrap, "iloc"):
        if columns is not None:
            data_to_wrap.columns = columns
        if index is not None:
            data_to_wrap.index = index
        return data_to_wrap

    pd = check_pandas_support("Setting output container to 'pandas'")
    return pd.DataFrame(data_to_wrap, index=index, columns=columns)


def get_output_config(method, estimator=None):
    """Get output configure based on estimator and global configuration.

    .. note:: Experimental API
        The `get_output_config` API is experimental and subject to change without
        deprecation.

    Parameters
    ----------
    method : {"transform"}
        Method to get container output for.

    estimator : estimator instance or None
        Estimator to get the output configuration from. If `None`, check global
        configuration is used.

    Returns
    -------
    config : dict
        Dictionary with keys:

        - "dense": specifies the dense container for `method`. This can be
          `"default"` or `"pandas"`.
    """
    est_sklearn_output_config = getattr(estimator, "_sklearn_output_config", {})
    if method in est_sklearn_output_config:
        dense_config = est_sklearn_output_config[method]
    else:
        dense_config = get_config()[f"{method}_output"]

    if dense_config not in {"default", "pandas"}:
        raise ValueError(
            f"output config must be 'default' or 'pandas' got {dense_config}"
        )

    return {"dense": dense_config}


def _wrap_data_with_container(method, data_to_wrap, original_input, estimator):
    """Wrap output with container based on an estimator's or global config.

    Parameters
    ----------
    method : {"transform"}
        Method to get container output for.

    data_to_wrap : ndarray
        Data to wrap with container.

    original_input : {ndarray, dataframe}
        Original input of function.

    index : array-like
        Index to attach to output.

    estimator : estimator instance
        Estimator to get the output configuration from.

    Returns
    -------
    wrapped_output : ndarray or DataFrame
        Wrapped output with column names and index or `output` itself if wrapping is
        not configured.
    """
    output_config = get_output_config(method, estimator)

    if output_config["dense"] == "default":
        return data_to_wrap

    # dense_config == "pandas"
    return _wrap_in_pandas_container(
        data_to_wrap=data_to_wrap,
        index=getattr(original_input, "index", None),
        columns=getattr(estimator, "get_feature_names_out", None),
    )


def _wrap_method_output(f, method):
    """Wrapper used by SetOutputMixin to automatically wrap methods."""

    @wraps(f)
    def wrapped(self, X, *args, **kwargs):
        data_to_wrap = f(self, X, *args, **kwargs)
        # only wrap the first output for cross decomposition
        if isinstance(data_to_wrap, tuple):
            return (
                _wrap_data_with_container(method, data_to_wrap[0], X, self),
                *data_to_wrap[1:],
            )

        return _wrap_data_with_container(method, data_to_wrap, X, self)

    return wrapped


def _auto_wrap_is_configured(self):
    return hasattr(self, "get_feature_names_out") and getattr(
        self, "_sklearn_auto_wrap_output", False
    )


class SetOutputMixin:
    """Mixin that dynamically wraps methods to return container based on config.

    Currently `SetOutputMixin` wraps `transform` and `fit_transform` and configures
    it based on `set_output` of the global configuration.

    `set_output` is only defined if `get_feature_names_out` is defined and
    `auto_wrap_output` is True.
    """

    def __init_subclass__(cls, auto_wrap_output=True, **kwargs):
        # Dynamically wraps `transform` and `fit_transform` and configure it's
        # output based on `set_output`.
        if not isinstance(auto_wrap_output, bool):
            raise ValueError("auto_wrap_output should be a bool or a callable")

        cls._sklearn_auto_wrap_output = auto_wrap_output
        if not auto_wrap_output:
            return

        if hasattr(cls, "transform"):
            cls.transform = _wrap_method_output(cls.transform, "transform")
        if hasattr(cls, "fit_transform"):
            cls.fit_transform = _wrap_method_output(cls.fit_transform, "transform")

    @available_if(_auto_wrap_is_configured)
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

    .. note:: Experimental API
        The `safe_set_output` API is experimental and subject to change without
        deprecation.

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
