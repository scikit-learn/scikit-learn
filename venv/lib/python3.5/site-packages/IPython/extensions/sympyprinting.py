"""
**DEPRECATED**

A print function that pretty prints sympy Basic objects.

:moduleauthor: Brian Granger

Usage
=====

Once the extension is loaded, Sympy Basic objects are automatically
pretty-printed.

As of SymPy 0.7.2, maintenance of this extension has moved to SymPy under
sympy.interactive.ipythonprinting, any modifications to account for changes to
SymPy should be submitted to SymPy rather than changed here. This module is
maintained here for backwards compatablitiy with old SymPy versions.

"""
#-----------------------------------------------------------------------------
#  Copyright (C) 2008  The IPython Development Team
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import warnings

def load_ipython_extension(ip):
    warnings.warn("The sympyprinting extension has moved to `sympy`, "
        "use `from sympy import init_printing; init_printing()`")
