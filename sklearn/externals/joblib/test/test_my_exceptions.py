"""
Test my automatically generate exceptions
"""
from nose.tools import assert_true

from .. import my_exceptions


def test_inheritance():
    assert_true(isinstance(my_exceptions.JoblibNameError(), NameError))
    assert_true(isinstance(my_exceptions.JoblibNameError(),
                            my_exceptions.JoblibException))
    assert_true(my_exceptions.JoblibNameError is
                my_exceptions._mk_exception(NameError)[0])
