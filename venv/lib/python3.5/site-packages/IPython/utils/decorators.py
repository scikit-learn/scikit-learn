# encoding: utf-8
"""Decorators that don't go anywhere else.

This module contains misc. decorators that don't really go with another module
in :mod:`IPython.utils`. Beore putting something here please see if it should
go into another topical module in :mod:`IPython.utils`.
"""

#-----------------------------------------------------------------------------
#  Copyright (C) 2008-2011  The IPython Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file COPYING, distributed as part of this software.
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------

def flag_calls(func):
    """Wrap a function to detect and flag when it gets called.

    This is a decorator which takes a function and wraps it in a function with
    a 'called' attribute. wrapper.called is initialized to False.

    The wrapper.called attribute is set to False right before each call to the
    wrapped function, so if the call fails it remains False.  After the call
    completes, wrapper.called is set to True and the output is returned.

    Testing for truth in wrapper.called allows you to determine if a call to
    func() was attempted and succeeded."""
    
    # don't wrap twice
    if hasattr(func, 'called'):
        return func

    def wrapper(*args,**kw):
        wrapper.called = False
        out = func(*args,**kw)
        wrapper.called = True
        return out

    wrapper.called = False
    wrapper.__doc__ = func.__doc__
    return wrapper

def undoc(func):
    """Mark a function or class as undocumented.
    
    This is found by inspecting the AST, so for now it must be used directly
    as @undoc, not as e.g. @decorators.undoc
    """
    return func

