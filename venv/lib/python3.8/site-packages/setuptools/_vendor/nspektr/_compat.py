import contextlib
import sys


if sys.version_info >= (3, 10):
    import importlib.metadata as metadata
else:
    import setuptools.extern.importlib_metadata as metadata  # type: ignore # noqa: F401


def repair_extras(extras):
    """
    Repair extras that appear as match objects.

    python/importlib_metadata#369 revealed a flaw in the EntryPoint
    implementation. This function wraps the extras to ensure
    they are proper strings even on older implementations.
    """
    with contextlib.suppress(AttributeError):
        return list(item.group(0) for item in extras)
    return extras
