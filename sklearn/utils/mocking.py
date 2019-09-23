from ._mocking import *  # noqa
from ..utils.deprecated import _raise_dep_warning_if_not_pytest

deprecated_path = 'sklearn.utils.mocking'
correct_path = 'sklearn.utils'

_raise_dep_warning_if_not_pytest(deprecated_path, correct_path)
