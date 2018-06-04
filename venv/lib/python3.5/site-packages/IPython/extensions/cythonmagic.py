# -*- coding: utf-8 -*-
"""
**DEPRECATED**

The cython magic has been integrated into Cython itself, 
which is now released in version 0.21.

cf github `Cython` organisation, `Cython` repo, under the 
file `Cython/Build/IpythonMagic.py`
"""
#-----------------------------------------------------------------------------
# Copyright (C) 2010-2011, IPython Development Team.
#-----------------------------------------------------------------------------

import warnings

## still load the magic in IPython 3.x, remove completely in future versions.
def load_ipython_extension(ip):
    """Load the extension in IPython."""

    warnings.warn("""The Cython magic has been moved to the Cython package""")
