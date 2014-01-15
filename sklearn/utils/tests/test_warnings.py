# Author: Olivier Grisel <olivier.grisel@ensta.org>
# License: BSD 3 clause

import warnings
from multiprocessing.pool import ThreadPool
from ..safe_warnings import catch_warnings


class MyCustomWarning(Warning):
    pass


class MyOtherCustomWarning(Warning):
    pass


def check_caught_warning(iteration=0):
    with catch_warnings(record_category=MyCustomWarning) as w:
        warnings.warn(MyCustomWarning())
    if len(w) != 1 or w[0].category != MyCustomWarning:
        raise AssertionError(
            "A MyCustomWarning should have been caught,"
            " got %r" % list(w))


def check_nested_caught_warnings(iteration=0):
    with catch_warnings(record_category=MyOtherCustomWarning) as w:
        check_caught_warning()
        warnings.warn(MyOtherCustomWarning())
    if len(w) != 1 or w[0].category != MyOtherCustomWarning:
        raise AssertionError(
            "A MyOtherCustomWarning should have been caught,"
            " got %r" % list(w))


def test_thread_safe_caught_warnings():
    # check the sequential behavior first for easier debugging
    check_caught_warning()
    check_nested_caught_warnings()

    # actual check for thread safety
    p = ThreadPool(16)
    p.map(check_caught_warning, range(5000))
    p.map(check_nested_caught_warnings, range(5000))
    p.close()
