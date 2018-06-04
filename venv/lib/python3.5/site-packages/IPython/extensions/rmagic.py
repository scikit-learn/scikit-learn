# -*- coding: utf-8 -*-

#-----------------------------------------------------------------------------
#  Copyright (C) 2012 The IPython Development Team
#-----------------------------------------------------------------------------

import warnings

def load_ipython_extension(ip):
    """Load the extension in IPython."""
    warnings.warn("The rmagic extension in IPython has moved to "
            "`rpy2.ipython`, please see `rpy2` documentation.")
