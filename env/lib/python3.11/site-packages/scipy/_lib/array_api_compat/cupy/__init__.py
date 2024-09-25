from cupy import * # noqa: F403

# from cupy import * doesn't overwrite these builtin names
from cupy import abs, max, min, round # noqa: F401

# These imports may overwrite names from the import * above.
from ._aliases import * # noqa: F403

# See the comment in the numpy __init__.py
__import__(__package__ + '.linalg')

__import__(__package__ + '.fft')

from ..common._helpers import * # noqa: F401,F403

__array_api_version__ = '2022.12'
