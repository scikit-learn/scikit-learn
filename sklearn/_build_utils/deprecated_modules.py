"""Generates submodule to allow deprecation of submodules and keeping git
blame."""
from pathlib import Path

# This is a set of 3-tuples consisting of
# (new_module_path, deprecated_path, correct_path)
_DEPRECATED_MODULES = {
    # TODO: Remove in 0.24
    ('_mocking', 'sklearn.utils.mocking', 'sklearn.utils')
}

_FILE_CONTENT_TEMPLATE = """from .{new_module_path} import *  # noqa
from {relative_dots}utils.deprecation import _raise_dep_warning_if_not_pytest

deprecated_path = '{deprecated_path}'
correct_path = '{correct_path}'

_raise_dep_warning_if_not_pytest(deprecated_path, correct_path)
"""


def _create_deprecated_modules_files():
    """Add submodules that will be deprecated. A file is created based
    on the deprecated submodule's name. When this submodule is imported a
    deprecation warning will be raised.
    """
    for new_module_path, deprecated_path, correct_path in _DEPRECATED_MODULES:
        relative_dots = deprecated_path.count(".") * "."
        deprecated_content = _FILE_CONTENT_TEMPLATE.format(
            new_module_path=new_module_path,
            relative_dots=relative_dots,
            deprecated_path=deprecated_path,
            correct_path=correct_path)
        deprecated_parts = deprecated_path.split(".")
        deprecated_parts[-1] = deprecated_parts[-1] + ".py"

        with Path(*deprecated_parts).open('w') as f:
            f.write(deprecated_content)
