from ._partial_dependence import *  # noqa
from ..utils.deprecation import _raise_dep_warning_if_not_pytest

# TODO: remove entire file in 0.24
deprecated_path = 'sklearn.ensemble.partial_dependence'
correct_path = 'sklearn.ensemble'

_raise_dep_warning_if_not_pytest(deprecated_path, correct_path)
