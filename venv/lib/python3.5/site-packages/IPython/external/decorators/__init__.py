try:
    from numpy.testing.decorators import *
    from numpy.testing.noseclasses import KnownFailure
except ImportError:
    from ._decorators import *
    try:
        from ._numpy_testing_noseclasses import KnownFailure
    except ImportError:
        pass
