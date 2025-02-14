"""Utilities for backend dispatching."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause


import functools
import os
import warnings
from importlib.metadata import entry_points


def _entry_points(group="sklearn_backends"):
    # Support Python versions before 3.10, which do not let you
    # filter groups directly.
    all_entry_points = entry_points()
    if hasattr(all_entry_points, "select"):
        selected_entry_points = all_entry_points.select(group=group)
    else:
        selected_entry_points = all_entry_points.get(group, ())
    return selected_entry_points


@functools.cache
def all_backends():
    """Find all backends that provide an implementation for an estimator"""
    backends = {}

    for backend in _entry_points():
        backends[backend.name] = backend.load()

    return backends


def dispatching_disabled():
    """Determine if dispatching has been disabled by the user"""
    no_dispatching = os.environ.get("SKLEARN_NO_DISPATCHING", False)
    if no_dispatching == "1":
        return True
    else:
        return False


def public_api_name(estimator):
    """Get the name of the public module for a scikit-learn estimator.

    This computes the name of the public submodule in which the estimator
    can be found.
    """
    full_name = estimator.__module__
    # This relies on the fact that scikit-learn does not use
    # sub-submodules in its public API.
    # This means that public name can be atmost `sklearn.foobar`.
    # XXX Without this assumption it is very hard (impossible?)
    # XXX to compute the public name without a lot of acrobatics
    # XXX of walking up/down namespaces.
    public_name = ".".join(full_name.split(".")[:2])
    return public_name


def find_backend(estimator, *fit_args, **fit_kwargs):
    """Find a suitable backend for the combination of estimator and arguments."""
    if dispatching_disabled():
        return

    # XXX Keeping it simple for now, we can re-use most of the infrastructure
    # XXX from the existing plugins PR to help find backends via entrypoints
    # XXX as well as deal with enabling this via a config option.
    # XXX For now using a hardcoded list of backends
    module_path = public_api_name(estimator)
    name = estimator.__class__.__name__
    estimator_name = f"{module_path}:{name}"

    # Backends are tried in alphabetical order, this makes things
    # predictable and stable across runs. Might need a better solution
    # when it becomes common that users have more than one backend
    # that would accept a particular call.
    for name in sorted(all_backends()):
        backend = all_backends()[name]

        wants_it = backend.can_has(estimator_name, estimator, *fit_args, **fit_kwargs)
        if not wants_it:
            continue

        Implementation = backend.get_implementation(estimator_name)

        impl = Implementation(estimator)
        warnings.warn(
            f"The '{estimator_name}' implementation was dispatched to"
            f" the '{name}' backend. Set SKLEARN_NO_DISPATCHING=1 to"
            " disable this behaviour.",
            DispatchNotification,
            # XXX from where should this warning originate?
            # XXX what is the right stacklevel here? 4 seems right for KMeans
            # XXX but in general?
            stacklevel=4,
        )
        return impl


class DispatchNotification(RuntimeWarning):
    """Notification issued when a estimator is dispatched to a backend."""

    pass
