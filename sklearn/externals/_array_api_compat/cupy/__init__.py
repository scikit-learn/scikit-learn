from cupy import *

# from cupy import * doesn't overwrite these builtin names
from cupy import abs, max, min, round

# These imports may overwrite names from the import * above.
from ._aliases import *

# See the comment in the numpy __init__.py
__import__(__package__ + '.linalg')

from .linalg import matrix_transpose, vecdot

from ..common._helpers import *

__array_api_version__ = '2021.12'
