import weakref
import gc

import pytest

from sklearn.utils._plot import _check_axes_has_been_used
from sklearn.utils._plot import _check_axes_has_display_object
from sklearn.utils._plot import _SKLEARN_AX_DISP_OBJ_REF_KEY


def test_axes_has_been_used(pyplot):
    fig, ax = pyplot.subplots()

    msg = "The ax was already used in a matplotlib plot function"
    _check_axes_has_been_used(ax)  # no error

    ax.plot([0, 1, 2], [1, 2, 3])
    with pytest.raises(ValueError, match=msg):
        _check_axes_has_been_used(ax)

    ax.clear()
    _check_axes_has_been_used(ax)  # no error

    ax.imshow([[0, 1, 2]])
    with pytest.raises(ValueError, match=msg):
        _check_axes_has_been_used(ax)


def test_check_axes_has_display_object(pyplot):
    fig, ax = pyplot.subplots()

    class MockDisplay1:
        pass

    display_1 = MockDisplay1()

    # first call sets display weakref in ax
    assert _check_axes_has_display_object(display_1, ax) == display_1
    assert hasattr(ax, _SKLEARN_AX_DISP_OBJ_REF_KEY)
    display_1_ref = getattr(ax, _SKLEARN_AX_DISP_OBJ_REF_KEY)
    assert isinstance(display_1_ref, weakref.ref)
    assert display_1_ref() == display_1

    class MockDisplay2:
        pass

    display_2 = MockDisplay2()

    # errors because ax already has a reference to a display object
    msg = ("The ax was already used by another display object which is not "
           "an instance of MockDisplay2")

    with pytest.raises(ValueError, match=msg):
        _check_axes_has_display_object(display_2, ax)

    # deleting the display_1 object
    del display_1
    gc.collect()  # python 3.5 needs to gc before the weak reference is removed
    assert display_1_ref() is None

    # With the display_1 deleted, a new display object can be added to ax
    assert _check_axes_has_display_object(display_2, ax) == display_2
    display_2_ref = getattr(ax, _SKLEARN_AX_DISP_OBJ_REF_KEY)
    assert isinstance(display_2_ref, weakref.ref)
    assert display_2_ref() == display_2
