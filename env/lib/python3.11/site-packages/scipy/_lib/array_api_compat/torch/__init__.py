from torch import * # noqa: F403

# Several names are not included in the above import *
import torch
for n in dir(torch):
    if (n.startswith('_')
        or n.endswith('_')
        or 'cuda' in n
        or 'cpu' in n
        or 'backward' in n):
        continue
    exec(n + ' = torch.' + n)

# These imports may overwrite names from the import * above.
from ._aliases import * # noqa: F403

# See the comment in the numpy __init__.py
__import__(__package__ + '.linalg')

__import__(__package__ + '.fft')

from ..common._helpers import * # noqa: F403

__array_api_version__ = '2022.12'
