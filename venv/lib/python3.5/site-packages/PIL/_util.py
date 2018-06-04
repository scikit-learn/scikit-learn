import os

if bytes is str:
    def isStringType(t):
        return isinstance(t, basestring)

    def isPath(f):
        return isinstance(f, basestring)
else:
    def isStringType(t):
        return isinstance(t, str)

    def isPath(f):
        return isinstance(f, (bytes, str))


# Checks if an object is a string, and that it points to a directory.
def isDirectory(f):
    return isPath(f) and os.path.isdir(f)


class deferred_error(object):
    def __init__(self, ex):
        self.ex = ex

    def __getattr__(self, elt):
        raise self.ex
