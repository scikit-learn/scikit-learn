"""
Various helper functions which are not part of the spec.

Functions which start with an underscore are for internal use only but helpers
that are in __all__ are intended as additional helper functions for use by end
users of the compat library.
"""
from __future__ import annotations

import sys
import math

def _is_numpy_array(x):
    # Avoid importing NumPy if it isn't already
    if 'numpy' not in sys.modules:
        return False

    import numpy as np

    # TODO: Should we reject ndarray subclasses?
    return isinstance(x, (np.ndarray, np.generic))

def _is_cupy_array(x):
    # Avoid importing NumPy if it isn't already
    if 'cupy' not in sys.modules:
        return False

    import cupy as cp

    # TODO: Should we reject ndarray subclasses?
    return isinstance(x, (cp.ndarray, cp.generic))

def _is_torch_array(x):
    # Avoid importing torch if it isn't already
    if 'torch' not in sys.modules:
        return False

    import torch

    # TODO: Should we reject ndarray subclasses?
    return isinstance(x, torch.Tensor)

def is_array_api_obj(x):
    """
    Check if x is an array API compatible array object.
    """
    return _is_numpy_array(x) \
        or _is_cupy_array(x) \
        or _is_torch_array(x) \
        or hasattr(x, '__array_namespace__')

def _check_api_version(api_version):
    if api_version is not None and api_version != '2021.12':
        raise ValueError("Only the 2021.12 version of the array API specification is currently supported")

def array_namespace(*xs, api_version=None, _use_compat=True):
    """
    Get the array API compatible namespace for the arrays `xs`.

    `xs` should contain one or more arrays.

    Typical usage is

        def your_function(x, y):
            xp = array_api_compat.array_namespace(x, y)
            # Now use xp as the array library namespace
            return xp.mean(x, axis=0) + 2*xp.std(y, axis=0)

    api_version should be the newest version of the spec that you need support
    for (currently the compat library wrapped APIs only support v2021.12).
    """
    namespaces = set()
    for x in xs:
        if isinstance(x, (tuple, list)):
            namespaces.add(array_namespace(*x, _use_compat=_use_compat))
        elif hasattr(x, '__array_namespace__'):
            namespaces.add(x.__array_namespace__(api_version=api_version))
        elif _is_numpy_array(x):
            _check_api_version(api_version)
            if _use_compat:
                from .. import numpy as numpy_namespace
                namespaces.add(numpy_namespace)
            else:
                import numpy as np
                namespaces.add(np)
        elif _is_cupy_array(x):
            _check_api_version(api_version)
            if _use_compat:
                from .. import cupy as cupy_namespace
                namespaces.add(cupy_namespace)
            else:
                import cupy as cp
                namespaces.add(cp)
        elif _is_torch_array(x):
            _check_api_version(api_version)
            if _use_compat:
                from .. import torch as torch_namespace
                namespaces.add(torch_namespace)
            else:
                import torch
                namespaces.add(torch)
        else:
            # TODO: Support Python scalars?
            raise TypeError("The input is not a supported array type")

    if not namespaces:
        raise TypeError("Unrecognized array input")

    if len(namespaces) != 1:
        raise TypeError(f"Multiple namespaces for array inputs: {namespaces}")

    xp, = namespaces

    return xp

# backwards compatibility alias
get_namespace = array_namespace

def _check_device(xp, device):
    if xp == sys.modules.get('numpy'):
        if device not in ["cpu", None]:
            raise ValueError(f"Unsupported device for NumPy: {device!r}")

# device() is not on numpy.ndarray and and to_device() is not on numpy.ndarray
# or cupy.ndarray. They are not included in array objects of this library
# because this library just reuses the respective ndarray classes without
# wrapping or subclassing them. These helper functions can be used instead of
# the wrapper functions for libraries that need to support both NumPy/CuPy and
# other libraries that use devices.
def device(x: "Array", /) -> "Device":
    """
    Hardware device the array data resides on.

    Parameters
    ----------
    x: array
        array instance from NumPy or an array API compatible library.

    Returns
    -------
    out: device
        a ``device`` object (see the "Device Support" section of the array API specification).
    """
    if _is_numpy_array(x):
        return "cpu"
    return x.device

# Based on cupy.array_api.Array.to_device
def _cupy_to_device(x, device, /, stream=None):
    import cupy as cp
    from cupy.cuda import Device as _Device
    from cupy.cuda import stream as stream_module
    from cupy_backends.cuda.api import runtime

    if device == x.device:
        return x
    elif not isinstance(device, _Device):
        raise ValueError(f"Unsupported device {device!r}")
    else:
        # see cupy/cupy#5985 for the reason how we handle device/stream here
        prev_device = runtime.getDevice()
        prev_stream: stream_module.Stream = None
        if stream is not None:
            prev_stream = stream_module.get_current_stream()
            # stream can be an int as specified in __dlpack__, or a CuPy stream
            if isinstance(stream, int):
                stream = cp.cuda.ExternalStream(stream)
            elif isinstance(stream, cp.cuda.Stream):
                pass
            else:
                raise ValueError('the input stream is not recognized')
            stream.use()
        try:
            runtime.setDevice(device.id)
            arr = x.copy()
        finally:
            runtime.setDevice(prev_device)
            if stream is not None:
                prev_stream.use()
        return arr

def _torch_to_device(x, device, /, stream=None):
    if stream is not None:
        raise NotImplementedError
    return x.to(device)

def to_device(x: "Array", device: "Device", /, *, stream: "Optional[Union[int, Any]]" = None) -> "Array":
    """
    Copy the array from the device on which it currently resides to the specified ``device``.

    Parameters
    ----------
    x: array
        array instance from NumPy or an array API compatible library.
    device: device
        a ``device`` object (see the "Device Support" section of the array API specification).
    stream: Optional[Union[int, Any]]
        stream object to use during copy. In addition to the types supported in ``array.__dlpack__``, implementations may choose to support any library-specific stream object with the caveat that any code using such an object would not be portable.

    Returns
    -------
    out: array
        an array with the same data and data type as ``x`` and located on the specified ``device``.

    .. note::
       If ``stream`` is given, the copy operation should be enqueued on the provided ``stream``; otherwise, the copy operation should be enqueued on the default stream/queue. Whether the copy is performed synchronously or asynchronously is implementation-dependent. Accordingly, if synchronization is required to guarantee data safety, this must be clearly explained in a conforming library's documentation.
    """
    if _is_numpy_array(x):
        if stream is not None:
            raise ValueError("The stream argument to to_device() is not supported")
        if device == 'cpu':
            return x
        raise ValueError(f"Unsupported device {device!r}")
    elif _is_cupy_array(x):
        # cupy does not yet have to_device
        return _cupy_to_device(x, device, stream=stream)
    elif _is_torch_array(x):
        return _torch_to_device(x, device, stream=stream)
    return x.to_device(device, stream=stream)

def size(x):
    """
    Return the total number of elements of x
    """
    if None in x.shape:
        return None
    return math.prod(x.shape)

__all__ = ['is_array_api_obj', 'array_namespace', 'get_namespace', 'device', 'to_device', 'size']
