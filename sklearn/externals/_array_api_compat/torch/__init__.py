from torch import *

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
from ._aliases import *

# See the comment in the numpy __init__.py
__import__(__package__ + '.linalg')

from ..common._helpers import *

__array_api_version__ = '2021.12'
