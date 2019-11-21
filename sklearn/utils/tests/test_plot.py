import pytest

from sklearn.utils._plot import _check_axes_has_been_used


def test_axes_has_been_used(pyplot):
    fig, ax = pyplot.subplots()

    msg = "The ax was already used in another plot function"
    _check_axes_has_been_used(ax)  # no error

    ax.plot([0, 1, 2], [1, 2, 3])
    with pytest.raises(ValueError, match=msg):
        _check_axes_has_been_used(ax)

    ax.clear()

    _check_axes_has_been_used(ax)  # no error

    ax.imshow([[0, 1, 2]])
    with pytest.raises(ValueError, match=msg):
        _check_axes_has_been_used(ax)
