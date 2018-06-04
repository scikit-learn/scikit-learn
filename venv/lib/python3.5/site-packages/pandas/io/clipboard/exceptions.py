import ctypes


class PyperclipException(RuntimeError):
    pass


class PyperclipWindowsException(PyperclipException):

    def __init__(self, message):
        message += " ({err})".format(err=ctypes.WinError())
        super(PyperclipWindowsException, self).__init__(message)
