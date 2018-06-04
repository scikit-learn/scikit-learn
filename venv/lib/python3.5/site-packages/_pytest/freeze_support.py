"""
Provides a function to report all internal modules for using freezing tools
pytest
"""
from __future__ import absolute_import, division, print_function


def freeze_includes():
    """
    Returns a list of module names used by py.test that should be
    included by cx_freeze.
    """
    import py
    import _pytest
    result = list(_iter_all_modules(py))
    result += list(_iter_all_modules(_pytest))
    return result


def _iter_all_modules(package, prefix=''):
    """
    Iterates over the names of all modules that can be found in the given
    package, recursively.
    Example:
        _iter_all_modules(_pytest) ->
            ['_pytest.assertion.newinterpret',
             '_pytest.capture',
             '_pytest.core',
             ...
            ]
    """
    import os
    import pkgutil
    if type(package) is not str:
        path, prefix = package.__path__[0], package.__name__ + '.'
    else:
        path = package
    for _, name, is_package in pkgutil.iter_modules([path]):
        if is_package:
            for m in _iter_all_modules(os.path.join(path, name), prefix=name + '.'):
                yield prefix + m
        else:
            yield prefix + name
