import sys
import warnings

import numpy as np
import pytest

from skimage._shared import testing
from skimage._shared.utils import (
    _supported_float_type,
    _validate_interpolation_order,
    change_default_value,
    channel_as_last_axis,
    check_nD,
    deprecate_func,
    deprecate_parameter,
    DEPRECATED,
)

complex_dtypes = [np.complex64, np.complex128]
if hasattr(np, 'complex256'):
    complex_dtypes += [np.complex256]

have_numpydoc = False
try:
    import numpydoc  # noqa: F401

    have_numpydoc = True
except ImportError:
    pass


def test_change_default_value():
    @change_default_value('arg1', new_value=-1, changed_version='0.12')
    def foo(arg0, arg1=0, arg2=1):
        """Expected docstring"""
        return arg0, arg1, arg2

    @change_default_value(
        'arg1',
        new_value=-1,
        changed_version='0.12',
        warning_msg="Custom warning message",
    )
    def bar(arg0, arg1=0, arg2=1):
        """Expected docstring"""
        return arg0, arg1, arg2

    # Assert warning messages
    with pytest.warns(FutureWarning) as record:
        assert foo(0) == (0, 0, 1)
        assert bar(0) == (0, 0, 1)

    expected_msg = (
        "The new recommended value for arg1 is -1. Until "
        "version 0.12, the default arg1 value is 0. From "
        "version 0.12, the arg1 default value will be -1. "
        "To avoid this warning, please explicitly set arg1 value."
    )

    assert str(record[0].message) == expected_msg
    assert str(record[1].message) == "Custom warning message"

    # Assert that nothing happens if arg1 is set
    with warnings.catch_warnings(record=True) as recorded:
        # No kwargs
        assert foo(0, 2) == (0, 2, 1)
        assert foo(0, arg1=0) == (0, 0, 1)

        # Function name and doc is preserved
        assert foo.__name__ == 'foo'
        if sys.flags.optimize < 2:
            # if PYTHONOPTIMIZE is set to 2, docstrings are stripped
            assert foo.__doc__ == 'Expected docstring'
    # Assert no warnings were raised
    assert len(recorded) == 0


def test_check_nD():
    z = np.random.random(200**2).reshape((200, 200))
    x = z[10:30, 30:10]
    with testing.raises(ValueError):
        check_nD(x, 2)


@pytest.mark.parametrize(
    'dtype', [bool, int, np.uint8, np.uint16, float, np.float32, np.float64]
)
@pytest.mark.parametrize('order', [None, -1, 0, 1, 2, 3, 4, 5, 6])
def test_validate_interpolation_order(dtype, order):
    if order is None:
        # Default order
        assert _validate_interpolation_order(dtype, None) == 0 if dtype == bool else 1
    elif order < 0 or order > 5:
        # Order not in valid range
        with testing.raises(ValueError):
            _validate_interpolation_order(dtype, order)
    elif dtype == bool and order != 0:
        # Deprecated order for bool array
        with pytest.raises(ValueError):
            _validate_interpolation_order(bool, order)
    else:
        # Valid use case
        assert _validate_interpolation_order(dtype, order) == order


@pytest.mark.parametrize(
    'dtype',
    [
        bool,
        np.float16,
        np.float32,
        np.float64,
        np.uint8,
        np.uint16,
        np.uint32,
        np.uint64,
        np.int8,
        np.int16,
        np.int32,
        np.int64,
    ],
)
def test_supported_float_dtype_real(dtype):
    float_dtype = _supported_float_type(dtype)
    if dtype in [np.float16, np.float32]:
        assert float_dtype == np.float32
    else:
        assert float_dtype == np.float64


@pytest.mark.parametrize('dtype', complex_dtypes)
@pytest.mark.parametrize('allow_complex', [False, True])
def test_supported_float_dtype_complex(dtype, allow_complex):
    if allow_complex:
        float_dtype = _supported_float_type(dtype, allow_complex=allow_complex)
        if dtype == np.complex64:
            assert float_dtype == np.complex64
        else:
            assert float_dtype == np.complex128
    else:
        with testing.raises(ValueError):
            _supported_float_type(dtype, allow_complex=allow_complex)


@pytest.mark.parametrize('dtype', ['f', 'float32', np.float32, np.dtype(np.float32)])
def test_supported_float_dtype_input_kinds(dtype):
    assert _supported_float_type(dtype) == np.float32


@pytest.mark.parametrize(
    'dtypes, expected',
    [
        ((np.float16, np.float64), np.float64),
        ((np.float32, np.uint16, np.int8), np.float64),
        ((np.float32, np.float16), np.float32),
    ],
)
def test_supported_float_dtype_sequence(dtypes, expected):
    float_dtype = _supported_float_type(dtypes)
    assert float_dtype == expected


@channel_as_last_axis(multichannel_output=False)
def _decorated_channel_axis_size(x, *, channel_axis=None):
    if channel_axis is None:
        return None
    assert channel_axis == -1
    return x.shape[-1]


@testing.parametrize('channel_axis', [None, 0, 1, 2, -1, -2, -3])
def test_decorated_channel_axis_shape(channel_axis):
    # Verify that channel_as_last_axis modifies the channel_axis as expected

    # need unique size per axis here
    x = np.zeros((2, 3, 4))

    size = _decorated_channel_axis_size(x, channel_axis=channel_axis)
    if channel_axis is None:
        assert size is None
    else:
        assert size == x.shape[channel_axis]


@deprecate_func(
    deprecated_version="x", removed_version="y", hint="You are on your own."
)
def _deprecated_func():
    """Dummy function used in `test_deprecate_func`.

    The decorated function must be outside the test function, otherwise it
    seems that the warning does not point at the calling location.
    """


def test_deprecate_func():
    with pytest.warns(FutureWarning) as record:
        _deprecated_func()
        testing.assert_stacklevel(record)

    assert len(record) == 1
    assert record[0].message.args[0] == (
        "`_deprecated_func` is deprecated since version x and will be removed in "
        "version y. You are on your own."
    )


@deprecate_parameter("old1", start_version="0.10", stop_version="0.12")
@deprecate_parameter("old0", start_version="0.10", stop_version="0.12")
def _func_deprecated_params(arg0, old0=DEPRECATED, old1=DEPRECATED, arg1=None):
    """Expected docstring.

    Parameters
    ----------
    arg0 : int
        First unchanged parameter.
    arg1 : int, optional
        Second unchanged parameter.
    """
    return arg0, old0, old1, arg1


@deprecate_parameter("old1", new_name="new0", start_version="0.10", stop_version="0.12")
@deprecate_parameter("old0", new_name="new1", start_version="0.10", stop_version="0.12")
def _func_replace_params(
    arg0, old0=DEPRECATED, old1=DEPRECATED, new0=None, new1=None, arg1=None
):
    """Expected docstring.

    Parameters
    ----------
    arg0 : int
        First unchanged parameter.
    new0 : int, optional
        First new parameter.

        .. versionadded:: 0.10
    new1 : int, optional
        Second new parameter.

        .. versionadded:: 0.10
    arg1 : int, optional
        Second unchanged parameter.
    """
    return arg0, old0, old1, new0, new1, arg1


class Test_deprecate_parameter:
    @pytest.mark.skipif(not have_numpydoc, reason="requires numpydoc")
    def test_docstring_removed_param(self):
        # function name and doc are preserved
        assert _func_deprecated_params.__name__ == "_func_deprecated_params"
        if sys.flags.optimize < 2:
            # if PYTHONOPTIMIZE is set to 2, docstrings are stripped
            assert (
                _func_deprecated_params.__doc__
                == """Expected docstring.


    Parameters
    ----------
    arg0 : int
        First unchanged parameter.
    arg1 : int, optional
        Second unchanged parameter.

    Other Parameters
    ----------------
    old0 : DEPRECATED
        `old0` is deprecated.

        .. deprecated:: 0.10
    old1 : DEPRECATED
        `old1` is deprecated.

        .. deprecated:: 0.10
"""
            )

    @pytest.mark.skipif(not have_numpydoc, reason="requires numpydoc")
    def test_docstring_replaced_param(self):
        assert _func_replace_params.__name__ == "_func_replace_params"
        if sys.flags.optimize < 2:
            # if PYTHONOPTIMIZE is set to 2, docstrings are stripped
            assert (
                _func_replace_params.__doc__
                == """Expected docstring.


    Parameters
    ----------
    arg0 : int
        First unchanged parameter.
    new0 : int, optional
        First new parameter.

        .. versionadded:: 0.10
    new1 : int, optional
        Second new parameter.

        .. versionadded:: 0.10
    arg1 : int, optional
        Second unchanged parameter.

    Other Parameters
    ----------------
    old0 : DEPRECATED
        Deprecated in favor of `new1`.

        .. deprecated:: 0.10
    old1 : DEPRECATED
        Deprecated in favor of `new0`.

        .. deprecated:: 0.10
"""
            )

    def test_warning_removed_param(self):
        match = (
            r".*`old[01]` is deprecated since version 0\.10 and will be removed "
            r"in 0\.12.* see the documentation of .*_func_deprecated_params`."
        )
        with pytest.warns(FutureWarning, match=match):
            assert _func_deprecated_params(1, 2) == (1, 2, DEPRECATED, None)
        with pytest.warns(FutureWarning, match=match):
            assert _func_deprecated_params(1, 2, 3) == (1, 2, 3, None)
        with pytest.warns(FutureWarning, match=match):
            assert _func_deprecated_params(1, old0=2) == (
                1,
                2,
                DEPRECATED,
                None,
            )
        with pytest.warns(FutureWarning, match=match):
            assert _func_deprecated_params(1, old1=2) == (
                1,
                DEPRECATED,
                2,
                None,
            )

        with warnings.catch_warnings(record=True) as record:
            assert _func_deprecated_params(1, arg1=3) == (1, DEPRECATED, DEPRECATED, 3)
        assert len(record) == 0

    def test_warning_replaced_param(self):
        match = (
            r".*`old[0,1]` is deprecated since version 0\.10 and will be removed "
            r"in 0\.12.* see the documentation of .*_func_replace_params`."
        )

        with pytest.warns(FutureWarning, match=match):
            assert _func_replace_params(1, 2) == (
                1,
                DEPRECATED,
                DEPRECATED,
                None,
                2,
                None,
            )

        with pytest.warns(FutureWarning, match=match) as records:
            assert _func_replace_params(1, 2, 3) == (
                1,
                DEPRECATED,
                DEPRECATED,
                3,
                2,
                None,
            )
        assert len(records) == 2
        assert "`old1` is deprecated" in records[0].message.args[0]
        assert "`old0` is deprecated" in records[1].message.args[0]

        with pytest.warns(FutureWarning, match=match):
            assert _func_replace_params(1, old0=2) == (
                1,
                DEPRECATED,
                DEPRECATED,
                None,
                2,
                None,
            )

        with pytest.warns(FutureWarning, match=match):
            assert _func_replace_params(1, old1=3) == (
                1,
                DEPRECATED,
                DEPRECATED,
                3,
                None,
                None,
            )

        # Otherwise, no warnings are emitted!
        with warnings.catch_warnings(record=True) as record:
            assert _func_replace_params(1, new0=2, new1=3) == (
                1,
                DEPRECATED,
                DEPRECATED,
                2,
                3,
                None,
            )
        assert len(record) == 0

    def test_missing_DEPRECATED(self):
        decorate = deprecate_parameter(
            "old", start_version="0.10", stop_version="0.12", stacklevel=2
        )

        def foo(arg0, old=None):
            return arg0, old

        with pytest.raises(RuntimeError, match="Expected .* <DEPRECATED>"):
            decorate(foo)

        def bar(arg0, old=DEPRECATED):
            return arg0

        assert decorate(bar)(1) == 1

    def test_new_keyword_only(self):
        @deprecate_parameter(
            "old",
            new_name="new",
            start_version="0.19",
            stop_version="0.21",
        )
        def foo(arg0, old=DEPRECATED, *, new=1, arg3=None):
            """Expected docstring"""
            return arg0, new, arg3

        # Assert that nothing happens when the function is called with the
        # new API
        with warnings.catch_warnings(record=True) as recorded:
            # No kwargs
            assert foo(0) == (0, 1, None)
            # Kwargs without deprecated argument
            assert foo(0, new=1, arg3=2) == (0, 1, 2)
            assert foo(0, new=2) == (0, 2, None)
            assert foo(0, arg3=2) == (0, 1, 2)
        assert len(recorded) == 0

    def test_conflicting_old_and_new(self):
        match = r".*`old[0,1]` is deprecated"
        with pytest.warns(FutureWarning, match=match):
            with pytest.raises(ValueError, match=".* avoid conflicting values"):
                _func_replace_params(1, old0=2, new1=2)

        with pytest.warns(FutureWarning, match=match):
            with pytest.raises(ValueError, match=".* avoid conflicting values"):
                _func_replace_params(1, old1=2, new0=2)

        with pytest.warns(FutureWarning, match=match):
            with pytest.raises(ValueError, match=".* avoid conflicting values"):
                _func_replace_params(1, old0=1, old1=1, new0=1, new1=1)

    def test_wrong_call_signature(self):
        """Check that normal errors for faulty calls are unchanged."""
        with pytest.raises(
            TypeError, match=r".* required positional argument\: 'arg0'"
        ):
            _func_replace_params()

        with pytest.warns(FutureWarning, match=r".*`old[0,1]` is deprecated"):
            with pytest.raises(
                TypeError, match=".* multiple values for argument 'old0'"
            ):
                _func_deprecated_params(1, 2, old0=2)

    def test_wrong_param_name(self):
        with pytest.raises(ValueError, match="'old' is not in list"):

            @deprecate_parameter("old", start_version="0.10", stop_version="0.12")
            def foo(arg0):
                pass

        with pytest.raises(ValueError, match="'new' is not in list"):

            @deprecate_parameter(
                "old", new_name="new", start_version="0.10", stop_version="0.12"
            )
            def bar(arg0, old, arg1):
                pass

    def test_warning_location(self):
        with pytest.warns(FutureWarning) as records:
            _func_deprecated_params(1, old0=2, old1=2)
            testing.assert_stacklevel(records)
        assert len(records) == 2

    def test_stacklevel(self):
        @deprecate_parameter(
            "old",
            start_version="0.19",
            stop_version="0.21",
        )
        def foo(arg0, old=DEPRECATED):
            pass

        with pytest.raises(RuntimeError, match="Set stacklevel manually"):
            foo(0, 1)

        @deprecate_parameter(
            "old",
            start_version="0.19",
            stop_version="0.21",
            stacklevel=2,
        )
        def bar(arg0, old=DEPRECATED):
            pass

        with pytest.warns(FutureWarning, match="`old` is deprecated") as records:
            bar(0, 1)
            testing.assert_stacklevel(records)
