from toolz.utils import raises


def test_raises():
    assert raises(ZeroDivisionError, lambda: 1 / 0)
    assert not raises(ZeroDivisionError, lambda: 1)
