"""Compiler tools with improved interactive support.

Provides compilation machinery similar to codeop, but with caching support so
we can provide interactive tracebacks.

Authors
-------
* Robert Kern
* Fernando Perez
* Thomas Kluyver
"""

# Note: though it might be more natural to name this module 'compiler', that
# name is in the stdlib and name collisions with the stdlib tend to produce
# weird problems (often with third-party tools).

#-----------------------------------------------------------------------------
#  Copyright (C) 2010-2011 The IPython Development Team.
#
#  Distributed under the terms of the BSD License.
#
#  The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

# Stdlib imports
import __future__
from ast import PyCF_ONLY_AST
import codeop
import functools
import hashlib
import linecache
import operator
import time

#-----------------------------------------------------------------------------
# Constants
#-----------------------------------------------------------------------------

# Roughtly equal to PyCF_MASK | PyCF_MASK_OBSOLETE as defined in pythonrun.h,
# this is used as a bitmask to extract future-related code flags.
PyCF_MASK = functools.reduce(operator.or_,
                             (getattr(__future__, fname).compiler_flag
                              for fname in __future__.all_feature_names))

#-----------------------------------------------------------------------------
# Local utilities
#-----------------------------------------------------------------------------

def code_name(code, number=0):
    """ Compute a (probably) unique name for code for caching.
    
    This now expects code to be unicode.
    """
    hash_digest = hashlib.sha1(code.encode("utf-8")).hexdigest()
    # Include the number and 12 characters of the hash in the name.  It's
    # pretty much impossible that in a single session we'll have collisions
    # even with truncated hashes, and the full one makes tracebacks too long
    return '<ipython-input-{0}-{1}>'.format(number, hash_digest[:12])

#-----------------------------------------------------------------------------
# Classes and functions
#-----------------------------------------------------------------------------

class CachingCompiler(codeop.Compile):
    """A compiler that caches code compiled from interactive statements.
    """

    def __init__(self):
        codeop.Compile.__init__(self)
        
        # This is ugly, but it must be done this way to allow multiple
        # simultaneous ipython instances to coexist.  Since Python itself
        # directly accesses the data structures in the linecache module, and
        # the cache therein is global, we must work with that data structure.
        # We must hold a reference to the original checkcache routine and call
        # that in our own check_cache() below, but the special IPython cache
        # must also be shared by all IPython instances.  If we were to hold
        # separate caches (one in each CachingCompiler instance), any call made
        # by Python itself to linecache.checkcache() would obliterate the
        # cached data from the other IPython instances.
        if not hasattr(linecache, '_ipython_cache'):
            linecache._ipython_cache = {}
        if not hasattr(linecache, '_checkcache_ori'):
            linecache._checkcache_ori = linecache.checkcache
        # Now, we must monkeypatch the linecache directly so that parts of the
        # stdlib that call it outside our control go through our codepath
        # (otherwise we'd lose our tracebacks).
        linecache.checkcache = check_linecache_ipython
        
    def ast_parse(self, source, filename='<unknown>', symbol='exec'):
        """Parse code to an AST with the current compiler flags active.
        
        Arguments are exactly the same as ast.parse (in the standard library),
        and are passed to the built-in compile function."""
        return compile(source, filename, symbol, self.flags | PyCF_ONLY_AST, 1)
    
    def reset_compiler_flags(self):
        """Reset compiler flags to default state."""
        # This value is copied from codeop.Compile.__init__, so if that ever
        # changes, it will need to be updated.
        self.flags = codeop.PyCF_DONT_IMPLY_DEDENT

    @property
    def compiler_flags(self):
        """Flags currently active in the compilation process.
        """
        return self.flags
        
    def cache(self, code, number=0):
        """Make a name for a block of code, and cache the code.
        
        Parameters
        ----------
        code : str
          The Python source code to cache.
        number : int
          A number which forms part of the code's name. Used for the execution
          counter.
          
        Returns
        -------
        The name of the cached code (as a string). Pass this as the filename
        argument to compilation, so that tracebacks are correctly hooked up.
        """
        name = code_name(code, number)
        entry = (len(code), time.time(),
                 [line+'\n' for line in code.splitlines()], name)
        linecache.cache[name] = entry
        linecache._ipython_cache[name] = entry
        return name

def check_linecache_ipython(*args):
    """Call linecache.checkcache() safely protecting our cached values.
    """
    # First call the original checkcache as intended
    linecache._checkcache_ori(*args)
    # Then, update back the cache with our data, so that tracebacks related
    # to our compiled codes can be produced.
    linecache.cache.update(linecache._ipython_cache)
