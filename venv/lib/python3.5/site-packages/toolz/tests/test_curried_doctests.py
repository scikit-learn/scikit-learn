import doctest
import toolz


def test_doctests():
    toolz.__test__ = {}
    for name, func in vars(toolz).items():
        if isinstance(func, toolz.curry):
            toolz.__test__[name] = func.func
    assert doctest.testmod(toolz).failed == 0
    del toolz.__test__
