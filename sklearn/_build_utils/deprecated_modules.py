"""Generates submodule to allow deprecation of submodules and keeping git
blame."""
from pathlib import Path


_DEPRECATED_MODULES = {
    # TODO: Remove in 0.24
    ('_mocking', 'sklearn.utils.mocking', 'sklearn.utils')
}

_DEPRECATE_TEMPLATE = """from .{module} import *  # noqa
from ..utils.deprecation import _raise_dep_warning_if_not_pytest

deprecated_path = '{deprecated_path}'
correct_path = '{correct_path}'

_raise_dep_warning_if_not_pytest(deprecated_path, correct_path)
"""


def _add_deprecated_submodules():
    """Add submodules that will be deprecated. A file is created based
    on the deprecated submodule's name. When this submodule is imported a
    deprecation warning will be raised.
    """
    for module, deprecated_path, correct_path in _DEPRECATED_MODULES:
        deprecated_content = _DEPRECATE_TEMPLATE.format(
            module=module, deprecated_path=deprecated_path,
            correct_path=correct_path)
        deprecated_parts = deprecated_path.split(".")
        deprecated_parts[-1] = deprecated_parts[-1] + ".py"

        with Path(*deprecated_parts).open('w') as f:
            f.write(deprecated_content)
