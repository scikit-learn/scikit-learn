"""[DEPRECATED] Utilities for connecting to kernels

Moved to IPython.kernel.connect
"""

import warnings
warnings.warn("IPython.lib.kernel moved to IPython.kernel.connect in IPython 1.0," 
        " and will be removed in IPython 6.0.",
    DeprecationWarning
)

from ipykernel.connect import *

