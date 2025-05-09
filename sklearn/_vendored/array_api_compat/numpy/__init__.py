from numpy import * # noqa: F403

# from numpy import * doesn't overwrite these builtin names
from numpy import abs, max, min, round # noqa: F401

# These imports may overwrite names from the import * above.
from ._aliases import * # noqa: F403

# Don't know why, but we have to do an absolute import to import linalg. If we
# instead do
#
# from . import linalg
#
# It doesn't overwrite np.linalg from above. The import is generated
# dynamically so that the library can be vendored.
__import__(__package__ + '.linalg')

__import__(__package__ + '.fft')

from .linalg import matrix_transpose, vecdot # noqa: F401

from ..common._helpers import * # noqa: F403

try:
    # Used in asarray(). Not present in older versions.
    from numpy import _CopyMode # noqa: F401
except ImportError:
    pass

__array_api_version__ = '2024.12'
