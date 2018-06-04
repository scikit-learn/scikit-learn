from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from .. import _copy_metadata
from . import knownfail
from .exceptions import KnownFailureDidNotFailTest


def knownfailureif(fail_condition, msg=None, known_exception_class=None):
    # based on numpy.testing.dec.knownfailureif
    if msg is None:
        msg = 'Test known to fail'

    def known_fail_decorator(f):
        def failer(*args, **kwargs):
            try:
                # Always run the test (to generate images).
                result = f(*args, **kwargs)
            except Exception as err:
                if fail_condition:
                    if known_exception_class is not None:
                        if not isinstance(err, known_exception_class):
                            # This is not the expected exception
                            raise
                    knownfail(msg)
                else:
                    raise
            if fail_condition and fail_condition != 'indeterminate':
                raise KnownFailureDidNotFailTest(msg)
            return result
        return _copy_metadata(f, failer)
    return known_fail_decorator
