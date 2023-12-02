import importlib
from functools import wraps
from typing import Protocol, runtime_checkable

import numpy as np
from scipy.sparse import issparse

from .._config import get_config
from ._available_if import available_if


def check_library_installed(library):
    """Check library is installed."""
    try:
        return importlib.import_module(library)
    except ImportError as exc:
        raise ImportError(
            f"Setting output container to '{library}' requires {library} to be"
            " installed"
        ) from exc


def get_columns(columns):
    if callable(columns):
        try:
            return columns()
        except Exception:
            return None
    return columns


@runtime_checkable
class ContainerAdapterProtocol(Protocol):
    container_lib: str

    def create_container(self, X_output, X_original, columns):
        """Create container from `X_output` with additional metadata.

        Parameters
        ----------
        X_output : {ndarray, dataframe}
            Data to wrap.

        X_original : {ndarray, dataframe}
            Original input dataframe. This is used to extract the metadata that should
            be passed to `X_output`, e.g. pandas row index.

        columns : callable, ndarray, or None
            The column names or a callable that returns the column names. The
            callable is useful if the column names require some computation. If `None`,
            then no columns are passed to the container's constructor.

        Returns
        -------
        wrapped_output : container_type
            `X_output` wrapped into the container type.
        """

    def is_supported_container(self, X):
        """Return True if X is a supported container.

        Parameters
        ----------
        Xs: container
            Containers to be checked.

        Returns
        -------
        is_supported_container : bool
            True if X is a supported container.
        """

    def rename_columns(self, X, columns):
        """Rename columns in `X`.

        Parameters
        ----------
        X : container
            Container which columns is updated.

        columns : ndarray of str
            Columns to update the `X`'s columns with.

        Returns
        -------
        updated_container : container
            Container with new names.
        """

    def hstack(self, Xs):
        """Stack containers horizontally (column-wise).

        Parameters
        ----------
        Xs : list of containers
            List of containers to stack.

        Returns
        -------
        stacked_Xs : container
            Stacked containers.
        """


class PandasAdapter:
    container_lib = "pandas"

    def create_container(self, X_output, X_original, columns):
        pd = check_library_installed("pandas")
        columns = get_columns(columns)
        index = X_original.index if isinstance(X_original, pd.DataFrame) else None

        if isinstance(X_output, pd.DataFrame):
            if columns is not None:
                X_output.columns = columns
            return X_output

        return pd.DataFrame(X_output, index=index, columns=columns, copy=False)

    def is_supported_container(self, X):
        pd = check_library_installed("pandas")
        return isinstance(X, pd.DataFrame)

    def rename_columns(self, X, columns):
        return X.rename(columns=dict(zip(X.columns, columns)))

    def hstack(self, Xs):
        pd = check_library_installed("pandas")
        return pd.concat(Xs, axis=1)


class PolarsAdapter:
    container_lib = "polars"

    def create_container(self, X_output, X_original, columns):
        pl = check_library_installed("polars")
        columns = get_columns(columns)

        if isinstance(columns, np.ndarray):
            columns = columns.tolist()

        if isinstance(X_output, pl.DataFrame):
            if columns is not None:
                return self.rename_columns(X_output, columns)
            return X_output

        return pl.DataFrame(X_output, schema=columns, orient="row")

    def is_supported_container(self, X):
        pl = check_library_installed("polars")
        return isinstance(X, pl.DataFrame)

    def rename_columns(self, X, columns):
        return X.rename(dict(zip(X.columns, columns)))

    def hstack(self, Xs):
        pl = check_library_installed("polars")
        return pl.concat(Xs, how="horizontal")


class ContainerAdaptersManager:
    def __init__(self):
        self.adapters = {}

    @property
    def supported_outputs(self):
        return {"default"} | set(self.adapters)

    def register(self, adapter):
        self.adapters[adapter.container_lib] = adapter


ADAPTERS_MANAGER = ContainerAdaptersManager()
ADAPTERS_MANAGER.register(PandasAdapter())
ADAPTERS_MANAGER.register(PolarsAdapter())


def _get_container_adapter(method, estimator=None):
    """Get container adapter."""
    dense_config = _get_output_config(method, estimator)["dense"]
    try:
        return ADAPTERS_MANAGER.adapters[dense_config]
    except KeyError:
        return None


def _get_output_config(method, estimator=None):
    """Get output config based on estimator and global configuration.

    Parameters
    ----------
    method : {"transform"}
        Estimator's method for which the output container is looked up.

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

    supported_outputs = ADAPTERS_MANAGER.supported_outputs
    if dense_config not in supported_outputs:
        raise ValueError(
            f"output config must be in {sorted(supported_outputs)}, got {dense_config}"
        )

    return {"dense": dense_config}


def _wrap_data_with_container(method, data_to_wrap, original_input, estimator):
    """Wrap output with container based on an estimator's or global config.

    Parameters
    ----------
    method : {"transform"}
        Estimator's method to get container output for.

    data_to_wrap : {ndarray, dataframe}
        Data to wrap with container.

    original_input : {ndarray, dataframe}
        Original input of function.

    estimator : estimator instance
        Estimator with to get the output configuration from.

    Returns
    -------
    output : {ndarray, dataframe}
        If the output config is "default" or the estimator is not configured
        for wrapping return `data_to_wrap` unchanged.
        If the output config is "pandas", return `data_to_wrap` as a pandas
        DataFrame.
    """
    output_config = _get_output_config(method, estimator)

    if output_config["dense"] == "default" or not _auto_wrap_is_configured(estimator):
        return data_to_wrap

    dense_config = output_config["dense"]
    if issparse(data_to_wrap):
        raise ValueError(
            "The transformer outputs a scipy sparse matrix. "
            "Try to set the transformer output to a dense array or disable "
            f"{dense_config.capitalize()} output with set_output(transform='default')."
        )

    adapter = ADAPTERS_MANAGER.adapters[dense_config]
    return adapter.create_container(
        data_to_wrap,
        original_input,
        columns=estimator.get_feature_names_out,
    )


def _wrap_method_output(f, method):
    """Wrapper used by `_SetOutputMixin` to automatically wrap methods."""

    @wraps(f)
    def wrapped(self, X, *args, **kwargs):
        data_to_wrap = f(self, X, *args, **kwargs)
        if isinstance(data_to_wrap, tuple):
            # only wrap the first output for cross decomposition
            return_tuple = (
                _wrap_data_with_container(method, data_to_wrap[0], X, self),
                *data_to_wrap[1:],
            )
            # Support for namedtuples `_make` is a documented API for namedtuples:
            # https://docs.python.org/3/library/collections.html#collections.somenamedtuple._make
            if hasattr(type(data_to_wrap), "_make"):
                return type(data_to_wrap)._make(return_tuple)
            return return_tuple

        return _wrap_data_with_container(method, data_to_wrap, X, self)

    return wrapped


def _auto_wrap_is_configured(estimator):
    """Return True if estimator is configured for auto-wrapping the transform method.

    `_SetOutputMixin` sets `_sklearn_auto_wrap_output_keys` to `set()` if auto wrapping
    is manually disabled.
    """
    auto_wrap_output_keys = getattr(estimator, "_sklearn_auto_wrap_output_keys", set())
    return (
        hasattr(estimator, "get_feature_names_out")
        and "transform" in auto_wrap_output_keys
    )


class _SetOutputMixin:
    """Mixin that dynamically wraps methods to return container based on config.

    Currently `_SetOutputMixin` wraps `transform` and `fit_transform` and configures
    it based on `set_output` of the global configuration.

    `set_output` is only defined if `get_feature_names_out` is defined and
    `auto_wrap_output_keys` is the default value.
    """

    def __init_subclass__(cls, auto_wrap_output_keys=("transform",), **kwargs):
        super().__init_subclass__(**kwargs)

        # Dynamically wraps `transform` and `fit_transform` and configure it's
        # output based on `set_output`.
        if not (
            isinstance(auto_wrap_output_keys, tuple) or auto_wrap_output_keys is None
        ):
            raise ValueError("auto_wrap_output_keys must be None or a tuple of keys.")

        if auto_wrap_output_keys is None:
            cls._sklearn_auto_wrap_output_keys = set()
            return

        # Mapping from method to key in configurations
        method_to_key = {
            "transform": "transform",
            "fit_transform": "transform",
            "inverse_transform": "transform",
        }
        cls._sklearn_auto_wrap_output_keys = set()

        for method, key in method_to_key.items():
            if not hasattr(cls, method) or key not in auto_wrap_output_keys:
                continue
            cls._sklearn_auto_wrap_output_keys.add(key)

            # Only wrap methods defined by cls itself
            if method not in cls.__dict__:
                continue
            wrapped_method = _wrap_method_output(getattr(cls, method), key)
            setattr(cls, method, wrapped_method)

    @available_if(_auto_wrap_is_configured)
    def set_output(self, *, transform=None):
        """Set output container.

        See :ref:`sphx_glr_auto_examples_miscellaneous_plot_set_output.py`
        for an example on how to use the API.

        Parameters
        ----------
        transform : {"default", "pandas"}, default=None
            Configure output of `transform` and `fit_transform`.

            - `"default"`: Default output format of a transformer
            - `"pandas"`: DataFrame output
            - `"polars"`: Polars output
            - `None`: Transform configuration is unchanged

            .. versionadded:: 1.4
                `"polars"` option was added.

        Returns
        -------
        self : estimator instance
            Estimator instance.
        """
        if transform is None:
            return self

        if not hasattr(self, "_sklearn_output_config"):
            self._sklearn_output_config = {}

        self._sklearn_output_config["transform"] = transform
        return self


def _safe_set_output(estimator, *, transform=None):
    """Safely call estimator.set_output and error if it not available.

    This is used by meta-estimators to set the output for child estimators.

    Parameters
    ----------
    estimator : estimator instance
        Estimator instance.

    transform : {"default", "pandas"}, default=None
        Configure output of the following estimator's methods:

        - `"transform"`
        - `"fit_transform"`

        If `None`, this operation is a no-op.

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
            "is not available."
        )
    return estimator.set_output(transform=transform)
