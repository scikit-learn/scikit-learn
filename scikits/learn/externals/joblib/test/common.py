"""
Small utilities for testing.
"""
import nose

# A decorator to run tests only when numpy is available
try:
    import numpy as np
    def with_numpy(func):
        """ A decorator to skip tests requiring numpy.
        """
        return func

except ImportError:    
    def with_numpy(func):
        """ A decorator to skip tests requiring numpy.
        """
        def my_func():
            raise nose.SkipTest('Test requires numpy')
        return my_func
    np = None


