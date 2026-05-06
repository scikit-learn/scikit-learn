import importlib
import inspect
import os
import pkgutil
import textwrap

import pytest

import sklearn
from sklearn import __version__
from sklearn.utils._openmp_helpers import _openmp_parallelism_enabled


@pytest.mark.thread_unsafe  # import side-effects
def test_extension_type_module():
    """Check that Cython extension types have a correct ``__module__``.

    When a subpackage containing Cython extension types has a misconfigured
    ``meson.build`` (e.g. missing ``__init__.py`` in its Cython tree), Cython
    cannot detect the package hierarchy and sets ``__module__`` to just the
    submodule name (e.g. ``'_loss'``) instead of the fully qualified
    ``'sklearn._loss._loss'``. This breaks downstream tools like skops that
    rely on ``__module__`` for serialization.
    """
    sklearn_path = [os.path.dirname(sklearn.__file__)]
    failures = []
    for _, modname, ispkg in pkgutil.walk_packages(
        path=sklearn_path, prefix="sklearn.", onerror=lambda _: None
    ):
        # Packages are directories, not modules that can hold extension
        # types. ``tests``, ``externals`` (vendored third-party code) and
        # ``_build_utils`` (build-time helpers that import ``Cython``, which
        # is not installed in the wheel test environment) are out of scope
        # for this check.
        if (
            ispkg
            or ".tests." in modname
            or ".externals." in modname
            or "._build_utils." in modname
        ):
            continue
        mod = importlib.import_module(modname)
        mod_file = getattr(mod, "__file__", "") or ""
        # Only compiled extension modules can produce the misconfigured
        # ``__module__`` this test guards against. Pure-Python modules get
        # the correct ``__module__`` from the import system by construction.
        if not mod_file.endswith((".so", ".pyd")):
            continue
        for name, cls in inspect.getmembers(mod, inspect.isclass):
            try:
                cls_file = inspect.getfile(cls)
            except TypeError:  # pragma: no cover
                # Raised for built-in types (``object``, stdlib C types) that
                # have no source file — they were not defined in ``mod``.
                continue  # pragma: no cover
            # Skip classes imported into ``mod`` from elsewhere (e.g. numpy,
            # scipy, or another sklearn module). Only classes whose source
            # file *is* this extension's .so are candidates for the bug.
            if cls_file != mod_file:
                continue
            if cls.__module__ != modname:
                failures.append(  # pragma: no cover
                    f"{modname}.{name}.__module__ == {cls.__module__!r}, "
                    f"expected {modname!r}"
                )
    assert not failures, "Extension types with incorrect __module__:\n" + "\n".join(
        failures
    )


def test_openmp_parallelism_enabled():
    # Check that sklearn is built with OpenMP-based parallelism enabled.
    # This test can be skipped by setting the environment variable
    # ``SKLEARN_SKIP_OPENMP_TEST``.
    if os.getenv("SKLEARN_SKIP_OPENMP_TEST"):
        pytest.skip("test explicitly skipped (SKLEARN_SKIP_OPENMP_TEST)")

    base_url = "dev" if __version__.endswith(".dev0") else "stable"
    err_msg = textwrap.dedent(
        """
        This test fails because scikit-learn has been built without OpenMP.
        This is not recommended since some estimators will run in sequential
        mode instead of leveraging thread-based parallelism.

        You can find instructions to build scikit-learn with OpenMP at this
        address:

            https://scikit-learn.org/{}/developers/advanced_installation.html

        You can skip this test by setting the environment variable
        SKLEARN_SKIP_OPENMP_TEST to any value.
        """
    ).format(base_url)

    assert _openmp_parallelism_enabled(), err_msg
