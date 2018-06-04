class UnpackException(Exception):
    """Deprecated.  Use Exception instead to catch all exception during unpacking."""


class BufferFull(UnpackException):
    pass


class OutOfData(UnpackException):
    pass


class UnpackValueError(UnpackException, ValueError):
    """Deprecated.  Use ValueError instead."""


class ExtraData(UnpackValueError):
    def __init__(self, unpacked, extra):
        self.unpacked = unpacked
        self.extra = extra

    def __str__(self):
        return "unpack(b) received extra data."


class PackException(Exception):
    """Deprecated.  Use Exception instead to catch all exception during packing."""


class PackValueError(PackException, ValueError):
    """PackValueError is raised when type of input data is supported but it's value is unsupported.

    Deprecated.  Use ValueError instead.
    """


class PackOverflowError(PackValueError, OverflowError):
    """PackOverflowError is raised when integer value is out of range of msgpack support [-2**31, 2**32).

    Deprecated.  Use ValueError instead.
    """
