"""Generates submodule to allow deprecation of submodules and keeping git
blame."""
from pathlib import Path
from contextlib import suppress

# This is a set of 3-tuples consisting of
# (new_module_name, deprecated_path, correct_import_path)
_DEPRECATED_MODULES = {
    # TODO: Remove in 0.24
    ('_mocking', 'sklearn.utils.mocking', 'sklearn.utils'),
    ('_bagging', 'sklearn.ensemble.bagging', 'sklearn.ensemble'),
    ('_base', 'sklearn.ensemble.base', 'sklearn.ensemble'),
    ('_forest', 'sklearn.ensemble.forest', 'sklearn.ensemble'),
    ('_gb', 'sklearn.ensemble.gradient_boosting', 'sklearn.ensemble'),
    ('_iforest', 'sklearn.ensemble.iforest', 'sklearn.ensemble'),
    ('_voting', 'sklearn.ensemble.voting', 'sklearn.ensemble'),
    ('_weight_boosting', 'sklearn.ensemble.weight_boosting', 'sklearn.ensemble'),
    ('_tree_base', 'sklearn.tree.tree', 'sklearn.tree'),
    ('_export', 'sklearn.tree.export', 'sklearn.tree'),
    ('_tree', 'sklearn.tree.tree_base', 'sklearn.tree'),
}

_FILE_CONTENT_TEMPLATE = """from .{new_module_name} import *  # noqa
from {relative_dots}utils.deprecation import _raise_dep_warning_if_not_pytest

deprecated_path = '{deprecated_path}'
correct_import_path = '{correct_import_path}'

_raise_dep_warning_if_not_pytest(deprecated_path, correct_import_path)
"""


def _get_deprecated_path(deprecated_path):
    deprecated_parts = deprecated_path.split(".")
    deprecated_parts[-1] = deprecated_parts[-1] + ".py"
    return Path(*deprecated_parts)


def _create_deprecated_modules_files():
    """Add submodules that will be deprecated. A file is created based
    on the deprecated submodule's name. When this submodule is imported a
    deprecation warning will be raised.
    """
    for (new_module_name, deprecated_path,
         correct_import_path) in _DEPRECATED_MODULES:
        relative_dots = deprecated_path.count(".") * "."
        deprecated_content = _FILE_CONTENT_TEMPLATE.format(
            new_module_name=new_module_name,
            relative_dots=relative_dots,
            deprecated_path=deprecated_path,
            correct_import_path=correct_import_path)

        with _get_deprecated_path(deprecated_path).open('w') as f:
            f.write(deprecated_content)


def _clean_deprecated_modules_files():
    """Removes submodules created by _create_deprecated_modules_files."""
    for (_, deprecated_path, _) in _DEPRECATED_MODULES:
        with suppress(FileNotFoundError):
            _get_deprecated_path(deprecated_path).unlink()


if __name__ == "__main__":
    _clean_deprecated_modules_files()
