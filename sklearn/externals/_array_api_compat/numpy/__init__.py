from numpy import *

# from numpy import * doesn't overwrite these builtin names
from numpy import abs, max, min, round

# These imports may overwrite names from the import * above.
from ._aliases import *

# Don't know why, but we have to do an absolute import to import linalg. If we
# instead do
#
# from . import linalg
#
# It doesn't overwrite np.linalg from above. The import is generated
# dynamically so that the library can be vendored.
__import__(__package__ + '.linalg')

from .linalg import matrix_transpose, vecdot

from ..common._helpers import *

__array_api_version__ = '2021.12'
