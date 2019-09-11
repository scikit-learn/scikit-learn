from ._tree_base import *  # noqa
from ..utils.deprecation import _raise_dep_warning_if_not_pytest

# TODO: remove entire file in 0.24
deprecated_path = 'sklearn.tree.tree'
correct_path = 'sklearn.tree'

_raise_dep_warning_if_not_pytest(deprecated_path, correct_path)
