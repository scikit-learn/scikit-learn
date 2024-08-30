import os
import threading
from contextlib import contextmanager as contextmanager

# Default configuration values
_global_config = {
    "assume_finite": bool(os.environ.get("SKLEARN_ASSUME_FINITE", False)),
    "working_memory": int(os.environ.get("SKLEARN_WORKING_MEMORY", 1024)),
    "print_changed_only": True,
    "display": "diagram",
    "pairwise_dist_chunk_size": int(os.environ.get("SKLEARN_PAIRWISE_DIST_CHUNK_SIZE", 256)),
    "enable_cython_pairwise_dist": True,
    "array_api_dispatch": False,
    "transform_output": "default",
    "enable_metadata_routing": False,
    "skip_parameter_validation": False,
}

_threadlocal = threading.local()


def _get_threadlocal_config():
    """Get a thread-local mutable configuration. If the configuration
    does not exist, copy the default global configuration."""
    if not hasattr(_threadlocal, "global_config"):
        _threadlocal.global_config = _global_config.copy()
    return _threadlocal.global_config


def get_config():
    """Retrieve current values for configuration set by :func:`set_config`.

    Returns
    -------
    config : dict
        Keys are parameter names that can be passed to :func:`set_config`.
    """
    return _get_threadlocal_config().copy()


def _apply_config_updates(local_config, **kwargs):
    """Apply configuration updates from kwargs to the local configuration."""
    for key, value in kwargs.items():
        if value is not None:
            local_config[key] = value


def set_config(
    assume_finite=None,
    working_memory=None,
    print_changed_only=None,
    display=None,
    pairwise_dist_chunk_size=None,
    enable_cython_pairwise_dist=None,
    array_api_dispatch=None,
    transform_output=None,
    enable_metadata_routing=None,
    skip_parameter_validation=None,
):
    """Set global scikit-learn configuration."""
    local_config = _get_threadlocal_config()

    config_updates = {
        "assume_finite": assume_finite,
        "working_memory": working_memory,
        "print_changed_only": print_changed_only,
        "display": display,
        "pairwise_dist_chunk_size": pairwise_dist_chunk_size,
        "enable_cython_pairwise_dist": enable_cython_pairwise_dist,
        "array_api_dispatch": array_api_dispatch,
        "transform_output": transform_output,
        "enable_metadata_routing": enable_metadata_routing,
        "skip_parameter_validation": skip_parameter_validation,
    }

    _apply_config_updates(local_config, **config_updates)

    if array_api_dispatch is not None:
        from .utils._array_api import _check_array_api_dispatch
        _check_array_api_dispatch(array_api_dispatch)


@contextmanager
def config_context(
    *,
    assume_finite=None,
    working_memory=None,
    print_changed_only=None,
    display=None,
    pairwise_dist_chunk_size=None,
    enable_cython_pairwise_dist=None,
    array_api_dispatch=None,
    transform_output=None,
    enable_metadata_routing=None,
    skip_parameter_validation=None,
):
    """Context manager for global scikit-learn configuration."""
    old_config = get_config()
    set_config(
        assume_finite=assume_finite,
        working_memory=working_memory,
        print_changed_only=print_changed_only,
        display=display,
        pairwise_dist_chunk_size=pairwise_dist_chunk_size,
        enable_cython_pairwise_dist=enable_cython_pairwise_dist,
        array_api_dispatch=array_api_dispatch,
        transform_output=transform_output,
        enable_metadata_routing=enable_metadata_routing,
        skip_parameter_validation=skip_parameter_validation,
    )

    try:
        yield
    finally:
        set_config(**old_config)
