#! /usr/bin/env python
# Last Change: Sat Jul 21 03:00 PM 2007 J

"""Mini script useful for top level testing of the learn package."""

from numpy.testing import NumpyTest

def additional_tests():
    # XXX: does this guarantee that the package is the one in the dev trunk, and
    # not scikits.foo installed somewhere else ?
    import scikits.learn
    np = NumpyTest(scikits.learn)
    return np._test_suite_from_all_tests(np.package, level = 10, verbosity = 1)
