from fnmatch import fnmatch


class GlobbingFilter(object):
    '''Filter module names using a set of globs.

    Objects are matched against the exclude list first, then the include list.
    Anything that passes through without matching either, is excluded.
    '''

    def __init__(self, include=None, exclude=None):
        if include is None and exclude is None:
            include = ['*']
            exclude = []
        elif include is None:
            include = ['*']
        elif exclude is None:
            exclude = []

        self.include = include
        self.exclude = exclude

    def __call__(self, full_name=None):
        for pattern in self.exclude:
            if fnmatch(full_name, pattern):
                return False

        for pattern in self.include:
            if fnmatch(full_name, pattern):
                return True

        return False
